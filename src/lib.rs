//! Anshel–Anshel–Goldfeld (AAG) style key exchange built over a non-commutative
//! algebraic platform: invertible quaternions modulo a prime `p`.
//!
//! Two parties, Alice and Bob, each hold a private *secret word* formed by
//! multiplying together a random sequence of generators drawn from their own
//! public generator set. They exchange conjugated versions of each other's
//! generator sets and, using their own secret pattern, each independently
//! telescopes the received conjugates down to a shared commutator value.
//! Because of how conjugation interacts with the secret words, both parties
//! arrive at (mutually inverse) values that hash to the same symmetric
//! session key.
//!
//! # Examples
//!
//! ```ignore
//! use qncke::PublicParams;
//! use rand_chacha::ChaCha8Rng;
//! use rand::SeedableRng;
//! use qncke::{Quaternion, quaternion_to_bytes};
//! let mut rng = ChaCha8Rng::seed_from_u64(42);
//! let params = PublicParams::generate(1019, 4, 4, 6, &mut rng);
//! let result = run_key_exchange_with_seeds(&params, Some(1), Some(2));
//! assert!(result.keys_match);
//! ```

#![allow(non_snake_case)]

use rand::Rng;
use sha3::{Digest, Sha3_256};

/// An element of the quaternion ring `(Z/pZ)[i, j, k]`, i.e. `a + b*i + c*j + d*k`
/// with all coefficients reduced modulo `p`.
///
/// Quaternions over a finite field form a non-commutative ring under the
/// standard Hamilton multiplication rules (reduced mod `p`), which makes them
/// a convenient platform group for Anshel–Anshel–Goldfeld-style protocols:
/// multiplication is cheap, but the underlying word/conjugacy problems are
/// non-trivial.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Quaternion {
    /// Real (scalar) component.
    pub a: u64,
    /// Coefficient of `i`.
    pub b: u64,
    /// Coefficient of `j`.
    pub c: u64,
    /// Coefficient of `k`.
    pub d: u64,
    /// Modulus defining the underlying ring `Z/pZ`. Should be prime for
    /// inversion to be well-defined.
    pub p: u64,
}

impl Quaternion {
    /// Create a new `Quaternion`, reducing all components modulo `p`.
    ///
    /// Components are accepted as `i128` so that negative intermediate
    /// values (e.g. from subtraction or conjugation) can be passed directly
    /// and correctly reduced into the range `[0, p)` via `rem_euclid`.
    ///
    /// # Examples
    ///
    /// ```ignore
    /// use qncke::{Quaternion, quaternion_to_bytes};
    /// let q = Quaternion::new(-1, 5, 0, 0, 7);
    /// assert_eq!(q.a, 6); // -1 mod 7 == 6
    /// assert_eq!(q.b, 5);
    /// ```
    pub fn new(a: i128, b: i128, c: i128, d: i128, p: u64) -> Self {
        let p_i = p as i128;
        Quaternion {
            a: a.rem_euclid(p_i) as u64,
            b: b.rem_euclid(p_i) as u64,
            c: c.rem_euclid(p_i) as u64,
            d: d.rem_euclid(p_i) as u64,
            p,
        }
    }

    /// Retrieve the four components as a tuple `(a, b, c, d)`.
    ///
    /// # Examples
    ///
    /// ```ignore
    /// use qncke::{Quaternion, quaternion_to_bytes};
    /// let q = Quaternion::new(1, 2, 3, 4, 11);
    /// assert_eq!(q.to_tuple(), (1, 2, 3, 4));
    /// ```
    pub fn to_tuple(&self) -> (u64, u64, u64, u64) {
        (self.a, self.b, self.c, self.d)
    }

    /// Algebraic conjugate: `q* = a - b*i - c*j - d*k` (mod `p`).
    ///
    /// This is the classical quaternion conjugate (negate the vector part),
    /// not to be confused with [`Quaternion::conjugate_by`], which performs
    /// group-theoretic conjugation `w^-1 * q * w`.
    ///
    /// # Examples
    ///
    /// ```ignore
    /// use qncke::{Quaternion, quaternion_to_bytes};
    /// let q = Quaternion::new(1, 2, 3, 4, 11);
    /// let conj = q.conjugate();
    /// assert_eq!(conj.to_tuple(), (1, 9, 8, 7)); // -2,-3,-4 mod 11
    /// ```
    pub fn conjugate(&self) -> Self {
        Quaternion::new(
            self.a as i128,
            -(self.b as i128),
            -(self.c as i128),
            -(self.d as i128),
            self.p,
        )
    }

    /// Multiplicative norm: `N(q) = a^2 + b^2 + c^2 + d^2` (mod `p`).
    ///
    /// The norm is multiplicative (`N(q1 * q2) = N(q1) * N(q2)` mod `p`) and
    /// is used to compute the multiplicative inverse of `q`.
    ///
    /// # Examples
    ///
    /// ```ignore
    /// use qncke::{Quaternion, quaternion_to_bytes};
    /// let q = Quaternion::new(1, 2, 0, 0, 11);
    /// assert_eq!(q.norm(), 5); // 1^2 + 2^2 = 5
    /// ```
    pub fn norm(&self) -> u64 {
        let a = self.a as u128;
        let b = self.b as u128;
        let c = self.c as u128;
        let d = self.d as u128;
        let p = self.p as u128;
        ((a * a + b * b + c * c + d * d) % p) as u64
    }

    /// Check whether this quaternion is invertible modulo `p`.
    ///
    /// A quaternion is invertible exactly when its [`norm`](Quaternion::norm)
    /// is non-zero modulo `p`.
    ///
    /// # Examples
    ///
    /// ```ignore
    /// use qncke::{Quaternion, quaternion_to_bytes};
    /// let zero = Quaternion::new(0, 0, 0, 0, 11);
    /// assert!(!zero.is_invertible());
    /// ```
    pub fn is_invertible(&self) -> bool {
        self.norm() % self.p != 0
    }

    /// Compute the multiplicative modular inverse: `q^-1 = q* * N(q)^-1` (mod `p`).
    ///
    /// Requires `p` to be prime and `N(q) != 0` (mod `p`); the norm's
    /// modular inverse is computed via Fermat's little theorem
    /// (`N(q)^(p-2) mod p`).
    ///
    /// # Panics
    ///
    /// Panics if the quaternion has zero norm modulo `p` (i.e. it is not
    /// invertible).
    ///
    /// # Examples
    ///
    /// ```ignore
    /// use qncke::{Quaternion, quaternion_to_bytes};
    /// let q = Quaternion::new(1, 2, 0, 0, 11);
    /// let inv = q.inverse();
    /// assert_eq!(q * inv, Quaternion::identity(11));
    /// ```
    pub fn inverse(&self) -> Self {
        let n = self.norm();
        if n == 0 {
            panic!("Quaternion has zero norm mod p; not invertible");
        }
        let n_inv = pow_mod(n, self.p - 2, self.p);
        let conj = self.conjugate();
        Quaternion::new(
            (conj.a as i128) * (n_inv as i128),
            (conj.b as i128) * (n_inv as i128),
            (conj.c as i128) * (n_inv as i128),
            (conj.d as i128) * (n_inv as i128),
            self.p,
        )
    }

    /// Group-theoretic conjugation of `self` by `w`: computes `w^-1 * self * w`.
    ///
    /// This is the operation underlying the AAG protocol's exchange step,
    /// distinct from the algebraic [`conjugate`](Quaternion::conjugate)
    /// (vector-part negation).
    ///
    /// # Examples
    ///
    /// ```ignore
    /// use qncke::{Quaternion, quaternion_to_bytes};
    /// let q = Quaternion::new(1, 2, 3, 4, 11);
    /// let w = Quaternion::new(2, 0, 1, 0, 11);
    /// let conjugated = q.conjugate_by(&w);
    /// assert_eq!(conjugated, w.inverse() * q * w);
    /// ```
    pub fn conjugate_by(&self, w: &Quaternion) -> Self {
        w.inverse() * *self * *w
    }

    /// The multiplicative identity quaternion `1 + 0i + 0j + 0k` over `Z/pZ`.
    ///
    /// # Examples
    ///
    /// ```ignore
    /// use qncke::{Quaternion, quaternion_to_bytes};
    /// let id = Quaternion::identity(11);
    /// assert_eq!(id.to_tuple(), (1, 0, 0, 0));
    /// ```
    pub fn identity(p: u64) -> Self {
        Quaternion::new(1, 0, 0, 0, p)
    }
}

/// Compute `base^exp mod modulus` using fast (binary) exponentiation.
///
/// Internal helper used by [`Quaternion::inverse`] to invert the norm via
/// Fermat's little theorem.
fn pow_mod(base: u64, mut exp: u64, modulus: u64) -> u64 {
    if modulus == 1 {
        return 0;
    }
    let mut result = 1u128;
    let mut base_128 = (base % modulus) as u128;
    let mod_128 = modulus as u128;
    while exp > 0 {
        if exp % 2 == 1 {
            result = (result * base_128) % mod_128;
        }
        base_128 = (base_128 * base_128) % mod_128;
        exp /= 2;
    }
    result as u64
}

impl std::ops::Add for Quaternion {
    type Output = Self;

    /// Component-wise addition modulo `p`.
    ///
    /// # Panics
    ///
    /// Panics if `self` and `other` are defined over different moduli.
    fn add(self, other: Self) -> Self {
        assert_eq!(self.p, other.p, "Quaternions defined over different moduli");
        Quaternion::new(
            (self.a as i128) + (other.a as i128),
            (self.b as i128) + (other.b as i128),
            (self.c as i128) + (other.c as i128),
            (self.d as i128) + (other.d as i128),
            self.p,
        )
    }
}

impl std::ops::Sub for Quaternion {
    type Output = Self;

    /// Component-wise subtraction modulo `p`.
    ///
    /// # Panics
    ///
    /// Panics if `self` and `other` are defined over different moduli.
    fn sub(self, other: Self) -> Self {
        assert_eq!(self.p, other.p, "Quaternions defined over different moduli");
        Quaternion::new(
            (self.a as i128) - (other.a as i128),
            (self.b as i128) - (other.b as i128),
            (self.c as i128) - (other.c as i128),
            (self.d as i128) - (other.d as i128),
            self.p,
        )
    }
}

impl std::ops::Mul for Quaternion {
    type Output = Self;

    /// Hamilton quaternion multiplication, reduced modulo `p`.
    ///
    /// Note that this operation is **non-commutative**: in general
    /// `q1 * q2 != q2 * q1`. This non-commutativity is what makes the
    /// quaternion ring a suitable platform for AAG-style protocols.
    ///
    /// # Panics
    ///
    /// Panics if `self` and `other` are defined over different moduli.
    fn mul(self, other: Self) -> Self {
        assert_eq!(self.p, other.p, "Quaternions defined over different moduli");
        let a1 = self.a as i128;
        let b1 = self.b as i128;
        let c1 = self.c as i128;
        let d1 = self.d as i128;

        let a2 = other.a as i128;
        let b2 = other.b as i128;
        let c2 = other.c as i128;
        let d2 = other.d as i128;

        let a = a1 * a2 - b1 * b2 - c1 * c2 - d1 * d2;
        let b = a1 * b2 + b1 * a2 + c1 * d2 - d1 * c2;
        let c = a1 * c2 - b1 * d2 + c1 * a2 + d1 * b2;
        let d = a1 * d2 + b1 * c2 - c1 * b2 + d1 * a2;

        Quaternion::new(a, b, c, d, self.p)
    }
}

/// Compute the left-to-right product of a sequence of quaternions.
///
/// Because quaternion multiplication is non-commutative, the order of the
/// slice matters: `product(&[q0, q1, q2])` computes `(q0 * q1) * q2`, not
/// any other association or ordering.
///
/// # Errors
///
/// Returns `Err("Cannot take product of empty sequence")` if `quaternions`
/// is empty.
///
/// # Examples
///
/// ```ignore
/// use qncke::{Quaternion, quaternion_to_bytes};
/// let p = 11;
/// let q0 = Quaternion::new(1, 1, 0, 0, p);
/// let q1 = Quaternion::new(0, 1, 1, 0, p);
/// assert_eq!(product(&[q0, q1]).unwrap(), q0 * q1);
/// ```
pub fn product(quaternions: &[Quaternion]) -> Result<Quaternion, &'static str> {
    if quaternions.is_empty() {
        return Err("Cannot take product of empty sequence");
    }
    let mut out = quaternions[0];
    for q in &quaternions[1..] {
        out = out * *q;
    }
    Ok(out)
}

/// Draw a uniformly random *invertible* quaternion over `Z/pZ`.
///
/// Repeatedly samples random components until an invertible quaternion
/// (non-zero norm mod `p`) is found. For reasonably large prime `p` this
/// terminates quickly, since most quaternions are invertible.
///
/// # Examples
///
/// ```ignore
/// use qncke::{Quaternion, quaternion_to_bytes};
/// use rand_chacha::ChaCha8Rng;
/// use rand::SeedableRng;
///
/// let mut rng = ChaCha8Rng::seed_from_u64(0);
/// let q = random_invertible_quaternion(1019, &mut rng);
/// assert!(q.is_invertible());
/// ```
pub fn random_invertible_quaternion<R: Rng>(p: u64, rng: &mut R) -> Quaternion {
    loop {
        let a = rng.gen_range(0..p) as i128;
        let b = rng.gen_range(0..p) as i128;
        let c = rng.gen_range(0..p) as i128;
        let d = rng.gen_range(0..p) as i128;
        let q = Quaternion::new(a, b, c, d, p);
        if q.is_invertible() {
            return q;
        }
    }
}

/// Fixed-width big-endian serialization of a quaternion's four `u64`
/// components into a 32-byte array (`a || b || c || d`).
///
/// # Examples
///
/// ```ignore
/// use qncke::{Quaternion, quaternion_to_bytes};
/// let q = Quaternion::new(1, 0, 0, 0, 11);
/// let bytes = quaternion_to_bytes(&q);
/// assert_eq!(bytes.len(), 32);
/// assert_eq!(&bytes[0..8], &1u64.to_be_bytes());
/// ```
pub fn quaternion_to_bytes(q: &Quaternion) -> [u8; 32] {
    let mut bytes = [0u8; 32];
    bytes[0..8].copy_from_slice(&q.a.to_be_bytes());
    bytes[8..16].copy_from_slice(&q.b.to_be_bytes());
    bytes[16..24].copy_from_slice(&q.c.to_be_bytes());
    bytes[24..32].copy_from_slice(&q.d.to_be_bytes());
    bytes
}

/// Compute the SHA3-256 digest of `data`, returned as a lowercase hex string.
///
/// # Examples
///
/// ```ignore
/// use qncke::{Quaternion, quaternion_to_bytes};
/// let digest = sha3_256_hex(b"hello world");
/// assert_eq!(digest.len(), 64); // 32 bytes -> 64 hex chars
/// ```
pub fn sha3_256_hex(data: &[u8]) -> String {
    let mut hasher = Sha3_256::new();
    hasher.update(data);
    let result = hasher.finalize();
    result.iter().map(|b| format!("{:02x}", b)).collect()
}

/// Derive a fixed-length symmetric session key (as a hex string) from a
/// shared quaternion value, by hashing its serialized bytes with SHA3-256.
///
/// # Examples
///
/// ```ignore
/// use qncke::{Quaternion, quaternion_to_bytes};
/// let q = Quaternion::new(1, 2, 3, 4, 1019);
/// let key = derive_session_key(&q);
/// assert_eq!(key.len(), 64);
/// ```
pub fn derive_session_key(q: &Quaternion) -> String {
    sha3_256_hex(&quaternion_to_bytes(q))
}

/// Expand a byte slice into its individual bits (0 or 1), most-significant
/// bit first within each byte.
///
/// # Examples
///
/// ```ignore
/// use qncke::{Quaternion, quaternion_to_bytes};
/// let bits = bytes_to_bits(&[0b1010_0000]);
/// assert_eq!(bits, vec![1, 0, 1, 0, 0, 0, 0, 0]);
/// ```
pub fn bytes_to_bits(data: &[u8]) -> Vec<u8> {
    let mut bits = Vec::with_capacity(data.len() * 8);
    for &byte in data {
        for i in (0..8).rev() {
            bits.push((byte >> i) & 1);
        }
    }
    bits
}

/// Compute the Hamming distance (number of differing bits) between two
/// equal-length hex-encoded byte strings.
///
/// # Errors
///
/// Returns an error if either string is not valid hex, or if the decoded
/// byte sequences have different lengths.
///
/// # Examples
///
/// ```ignore
/// use qncke::{Quaternion, quaternion_to_bytes};
/// assert_eq!(hamming_distance_hex("00", "ff").unwrap(), 8);
/// assert_eq!(hamming_distance_hex("0f", "00").unwrap(), 4);
/// ```
pub fn hamming_distance_hex(h1: &str, h2: &str) -> Result<usize, &'static str> {
    let b1 = hex_decode(h1)?;
    let b2 = hex_decode(h2)?;
    if b1.len() != b2.len() {
        return Err("Byte sequences must be of equal length");
    }
    let mut dist = 0;
    for (x, y) in b1.iter().zip(b2.iter()) {
        dist += (x ^ y).count_ones() as usize;
    }
    Ok(dist)
}

/// Decode a hex string into raw bytes.
///
/// Internal helper used by [`hamming_distance_hex`].
///
/// # Errors
///
/// Returns an error if `s` has odd length or contains non-hex characters.
fn hex_decode(s: &str) -> Result<Vec<u8>, &'static str> {
    if s.len() % 2 != 0 {
        return Err("Hex string must have an even length");
    }
    let mut bytes = Vec::with_capacity(s.len() / 2);
    for i in (0..s.len()).step_by(2) {
        let byte_str = &s[i..i + 2];
        let byte = u8::from_str_radix(byte_str, 16).map_err(|_| "Invalid hex character")?;
        bytes.push(byte);
    }
    Ok(bytes)
}

/// Public parameters shared by both parties before running the key exchange.
///
/// Defines the modulus `p`, each party's public generator set (`A` for
/// Alice, `B` for Bob), and the secret-word length `L` used by both parties.
#[derive(Debug, Clone)]
pub struct PublicParams {
    /// Prime modulus defining the quaternion ring `Z/pZ`.
    pub p: u64,
    /// Alice's public generator set.
    pub A: Vec<Quaternion>,
    /// Bob's public generator set.
    pub B: Vec<Quaternion>,
    /// Length of the secret word (number of generators multiplied together)
    /// each party uses.
    pub L: usize,
}

impl PublicParams {
    /// Generate public parameters: a prime modulus `p`, `k` random
    /// invertible generators for Alice, `m` random invertible generators
    /// for Bob, and a shared secret-word length `L`.
    ///
    /// # Examples
    ///
    /// ```ignore
    /// use qncke::{Quaternion, quaternion_to_bytes};
    /// use rand_chacha::ChaCha8Rng;
    /// use rand::SeedableRng;
    ///
    /// let mut rng = ChaCha8Rng::seed_from_u64(7);
    /// let params = PublicParams::generate(1019, 4, 4, 6, &mut rng);
    /// assert_eq!(params.A.len(), 4);
    /// assert_eq!(params.B.len(), 4);
    /// ```
    pub fn generate<R: Rng>(p: u64, k: usize, m: usize, L: usize, rng: &mut R) -> Self {
        let mut A = Vec::with_capacity(k);
        for _ in 0..k {
            A.push(random_invertible_quaternion(p, rng));
        }
        let mut B = Vec::with_capacity(m);
        for _ in 0..m {
            B.push(random_invertible_quaternion(p, rng));
        }
        PublicParams { p, A, B, L }
    }
}

/// A single participant in the AAG key exchange protocol.
///
/// Holds a private random `pattern` (indices into `generators`) and the
/// resulting `secret_word` (the product of the chosen generators), in
/// addition to its own public generator set.
#[derive(Debug, Clone)]
pub struct Party {
    /// Human-readable name for this party (e.g. `"Alice"`, `"Bob"`).
    pub name: String,
    /// This party's own public generator set.
    pub generators: Vec<Quaternion>,
    /// Length of the secret word.
    pub L: usize,
    /// Private sequence of indices into `generators` defining the secret
    /// word. Must never be revealed to the other party.
    pub pattern: Vec<usize>,
    /// The party's private secret word: the left-to-right product of
    /// `generators[pattern[i]]` for each `i`.
    pub secret_word: Quaternion,
}

impl Party {
    /// Create a new `Party`, immediately generating a random private
    /// `pattern` of length `L` and the corresponding `secret_word`.
    ///
    /// # Panics
    ///
    /// Panics if `L == 0` (the empty product is undefined; see
    /// [`product`]).
    ///
    /// # Examples
    ///
    /// ```ignore
    /// use qncke::{Quaternion, quaternion_to_bytes};
    /// use rand_chacha::ChaCha8Rng;
    /// use rand::SeedableRng;
    ///
    /// let mut rng = ChaCha8Rng::seed_from_u64(1);
    /// let generators = vec![Quaternion::new(1, 1, 0, 0, 1019)];
    /// let alice = Party::new("Alice", generators, 4, &mut rng);
    /// assert_eq!(alice.pattern.len(), 4);
    /// ```
    pub fn new<R: Rng>(name: &str, generators: Vec<Quaternion>, L: usize, rng: &mut R) -> Self {
        let n_gen = generators.len();
        let mut pattern = Vec::with_capacity(L);
        for _ in 0..L {
            pattern.push(rng.gen_range(0..n_gen));
        }
        let secret_word_components: Vec<Quaternion> = pattern.iter().map(|&idx| generators[idx]).collect();
        let secret_word = product(&secret_word_components).expect("L > 0");
        Party {
            name: name.to_string(),
            generators,
            L,
            pattern,
            secret_word,
        }
    }

    /// Round-1 message: conjugate each element of the other party's
    /// generator set by this party's secret word, i.e. compute
    /// `w^-1 * g * w` for every `g` in `other_generators`.
    ///
    /// This is the value that gets sent publicly to the other party.
    ///
    /// # Examples
    ///
    /// ```ignore
    /// use qncke::{Quaternion, quaternion_to_bytes};
    /// use rand_chacha::ChaCha8Rng;
    /// use rand::SeedableRng;
    ///
    /// let mut rng = ChaCha8Rng::seed_from_u64(2);
    /// let gens = vec![Quaternion::new(1, 1, 0, 0, 1019)];
    /// let alice = Party::new("Alice", gens.clone(), 3, &mut rng);
    /// let bob_gens = vec![Quaternion::new(0, 1, 1, 0, 1019)];
    /// let conjugated = alice.conjugate_set(&bob_gens);
    /// assert_eq!(conjugated.len(), bob_gens.len());
    /// ```
    pub fn conjugate_set(&self, other_generators: &[Quaternion]) -> Vec<Quaternion> {
        let w_inv = self.secret_word.inverse();
        other_generators
            .iter()
            .map(|&g| w_inv * g * self.secret_word)
            .collect()
    }

    /// Round-2: telescope the received conjugated generators using this
    /// party's own private `pattern`, then left-multiply by `w^-1` to
    /// recover the AAG commutator value shared with the other party.
    ///
    /// `received_conjugates` must be indexable by every value in
    /// `self.pattern` (i.e. it should have the same length as the other
    /// party's *own* generator set, since that is what was conjugated).
    ///
    /// # Panics
    ///
    /// Panics (via out-of-bounds indexing) if `received_conjugates` is
    /// shorter than `self.pattern` requires.
    pub fn derive_shared(&self, received_conjugates: &[Quaternion]) -> Quaternion {
        let telescoped_components: Vec<Quaternion> = self.pattern.iter().map(|&idx| received_conjugates[idx]).collect();
        let telescoped = product(&telescoped_components).expect("L > 0");
        self.secret_word.inverse() * telescoped
    }
}

/// The full output of a completed (or attempted) AAG key exchange between
/// Alice and Bob.
#[derive(Debug, Clone)]
pub struct ExchangeResult {
    /// Alice's private pattern (revealed here for testing/debugging only;
    /// in a real protocol this would never leave Alice's side).
    pub alice_pattern: Vec<usize>,
    /// Bob's private pattern (likewise, for testing/debugging only).
    pub bob_pattern: Vec<usize>,
    /// The commutator value Alice derives.
    pub K_alice: Quaternion,
    /// The raw (pre-normalization) commutator value Bob derives.
    pub K_bob_raw: Quaternion,
    /// Bob's value after the public normalization step (inversion), which
    /// should equal `K_alice` on success.
    pub K_bob_normalized: Quaternion,
    /// Whether `K_alice == K_bob_normalized`, i.e. whether the exchange
    /// succeeded.
    pub keys_match: bool,
    /// The derived session key (hex-encoded SHA3-256 digest), present only
    /// if `keys_match` is `true`.
    pub session_key_hex: Option<String>,
}

/// Run a full two-party AAG key exchange using the given public parameters
/// and per-party random number generators.
///
/// Generates fresh `Party` instances for Alice and Bob (each with its own
/// random secret pattern and secret word), performs the round-1 conjugated
/// generator-set exchange, and has each party independently derive the
/// shared commutator value in round 2. Returns the full [`ExchangeResult`],
/// including whether the keys matched and, if so, the resulting session key.
///
/// # Examples
///
/// ```ignore
/// use qncke::{Quaternion, quaternion_to_bytes};
/// use rand_chacha::ChaCha8Rng;
/// use rand::SeedableRng;
///
/// let mut rng_params = ChaCha8Rng::seed_from_u64(123);
/// let params = PublicParams::generate(1019, 4, 4, 6, &mut rng_params);
///
/// let mut rng_alice = ChaCha8Rng::seed_from_u64(1);
/// let mut rng_bob = ChaCha8Rng::seed_from_u64(2);
/// let result = run_key_exchange(&params, &mut rng_alice, &mut rng_bob);
/// assert!(result.keys_match);
/// ```
pub fn run_key_exchange<R1: Rng, R2: Rng>(
    params: &PublicParams,
    rng_alice: &mut R1,
    rng_bob: &mut R2,
) -> ExchangeResult {
    let alice = Party::new("Alice", params.A.clone(), params.L, rng_alice);
    let bob = Party::new("Bob", params.B.clone(), params.L, rng_bob);

    // Round 1 — public exchange of conjugated generator sets
    let c = alice.conjugate_set(&params.B); // Alice -> Bob
    let d = bob.conjugate_set(&params.A);   // Bob -> Alice

    // Round 2 — each side derives the (mutually inverse) commutator value
    let k_alice = alice.derive_shared(&d);
    let k_bob_raw = bob.derive_shared(&c);

    // Public normalization step: invert Bob's raw key
    let k_bob_normalized = k_bob_raw.inverse();

    let keys_match = k_alice == k_bob_normalized;
    let session_key_hex = if keys_match {
        Some(derive_session_key(&k_alice))
    } else {
        None
    };

    ExchangeResult {
        alice_pattern: alice.pattern,
        bob_pattern: bob.pattern,
        K_alice: k_alice,
        K_bob_raw: k_bob_raw,
        K_bob_normalized: k_bob_normalized,
        keys_match,
        session_key_hex,
    }
}

use rand::SeedableRng;
use rand_chacha::ChaCha8Rng;

/// Convenience wrapper around [`run_key_exchange`] that seeds each party's
/// RNG from an optional `u64` seed (falling back to OS entropy when `None`
/// is given), using the deterministic [`ChaCha8Rng`] generator.
///
/// Useful for reproducible tests and demonstrations, since identical seeds
/// always produce identical patterns, secret words, and (on success)
/// session keys.
///
/// # Examples
///
/// ```ignore
/// use qncke::{Quaternion, quaternion_to_bytes};
/// use rand_chacha::ChaCha8Rng;
/// use rand::SeedableRng;
/// let mut rng_params = rand_chacha::ChaCha8Rng::seed_from_u64(0);
/// use rand::SeedableRng;
/// let params = PublicParams::generate(1019, 4, 4, 6, &mut rng_params);
///
/// let result = run_key_exchange_with_seeds(&params, Some(10), Some(20));
/// assert!(result.keys_match);
/// ```
pub fn run_key_exchange_with_seeds(
    params: &PublicParams,
    seed_alice: Option<u64>,
    seed_bob: Option<u64>,
) -> ExchangeResult {
    let mut rng_alice = match seed_alice {
        Some(s) => ChaCha8Rng::seed_from_u64(s),
        None => ChaCha8Rng::from_entropy(),
    };
    let mut rng_bob = match seed_bob {
        Some(s) => ChaCha8Rng::seed_from_u64(s),
        None => ChaCha8Rng::from_entropy(),
    };
    run_key_exchange(params, &mut rng_alice, &mut rng_bob)
}

#[cfg(test)]
mod tests {
    use super::*;

    const P: u64 = 1019;

    #[test]
    fn test_associativity() {
        let mut rng = ChaCha8Rng::seed_from_u64(0);
        for _ in 0..20 {
            let a = random_invertible_quaternion(P, &mut rng);
            let b = random_invertible_quaternion(P, &mut rng);
            let c = random_invertible_quaternion(P, &mut rng);
            assert_eq!((a * b) * c, a * (b * c));
        }
    }

    #[test]
    fn test_inverse_correctness() {
        let mut rng = ChaCha8Rng::seed_from_u64(1);
        let identity = Quaternion::identity(P);
        for _ in 0..20 {
            let q = random_invertible_quaternion(P, &mut rng);
            assert_eq!(q * q.inverse(), identity);
            assert_eq!(q.inverse() * q, identity);
        }
    }

    #[test]
    fn test_norm_multiplicativity() {
        let mut rng = ChaCha8Rng::seed_from_u64(2);
        for _ in 0..20 {
            let a = random_invertible_quaternion(P, &mut rng);
            let b = random_invertible_quaternion(P, &mut rng);
            assert_eq!((a * b).norm(), (a.norm() * b.norm()) % P);
        }
    }

    #[test]
    fn test_product_helper_matches_manual_multiplication() {
        let mut rng = ChaCha8Rng::seed_from_u64(3);
        let qs = vec![
            random_invertible_quaternion(P, &mut rng),
            random_invertible_quaternion(P, &mut rng),
            random_invertible_quaternion(P, &mut rng),
            random_invertible_quaternion(P, &mut rng),
        ];
        let manual = qs[0] * qs[1] * qs[2] * qs[3];
        assert_eq!(product(&qs).unwrap(), manual);
    }

    #[test]
    fn test_key_agreement_success() {
        let mut rng_params = ChaCha8Rng::seed_from_u64(123);
        let params = PublicParams::generate(P, 4, 4, 6, &mut rng_params);
        for trial in 0..15 {
            let result = run_key_exchange_with_seeds(&params, Some(trial), Some(trial + 1000));
            assert!(result.keys_match, "Shared keys diverged on trial {}", trial);
            let key_hex = result.session_key_hex.expect("Session key should be present");
            assert_eq!(key_hex.len(), 64);
        }
    }

    #[test]
    fn test_key_agreement_changes_with_different_secrets() {
        let mut rng_params = ChaCha8Rng::seed_from_u64(321);
        let params = PublicParams::generate(P, 4, 4, 6, &mut rng_params);
        let r1 = run_key_exchange_with_seeds(&params, Some(1), Some(2));
        let r2 = run_key_exchange_with_seeds(&params, Some(3), Some(4));
        assert!(r1.keys_match);
        assert!(r2.keys_match);
        assert_ne!(r1.session_key_hex.unwrap(), r2.session_key_hex.unwrap());
    }
}