# QNCKE

A small Rust library implementing an **Anshel–Anshel–Goldfeld (AAG) style key exchange** over a non-commutative platform: **invertible quaternions modulo a prime `p`**.

It's a from-scratch, dependency-light playground for exploring non-commutative key exchange — group-theoretic conjugation, commutator derivation, and session-key hashing — built on a hand-rolled finite quaternion ring rather than braid groups or matrix groups (the more typical AAG platforms).

> ⚠️ **Educational / experimental code.** This is not a vetted cryptographic primitive — see [Security Status](#security-status) before using it for anything beyond learning and experimentation.

## How it works

Two parties, **Alice** and **Bob**, want to agree on a shared secret over a public channel.

1. **Setup.** Public parameters fix a prime modulus `p`, a generator set `A` (public, for Alice) and `B` (public, for Bob), and a secret-word length `L`.
2. **Private secret word.** Each party privately picks a random sequence (`pattern`) of `L` indices into *their own* generator set and multiplies the corresponding generators together (left to right) to get a private `secret_word` (`w`).
3. **Round 1 — conjugate and send.** Alice conjugates every element of Bob's public set `B` by her own secret word (`w⁻¹ · g · w` for each `g`), and sends the results to Bob. Bob does the same with Alice's set `A` and sends the results to Alice.
4. **Round 2 — telescope and derive.** Each party uses their own private `pattern` to pick out and multiply together the conjugated values they received, then left-multiplies by the inverse of their own secret word. This produces an AAG commutator-style value.
5. **Normalize and compare.** Bob's derived value is inverted as a public normalization step. If the protocol is implemented correctly, Alice's value and Bob's normalized value are equal — both sides now hold the same quaternion, which is hashed (SHA3-256) into a shared symmetric session key.

The quaternion ring is used here purely as a non-commutative algebraic stand-in — multiplication is cheap to compute, but order matters (`q1 * q2 ≠ q2 * q1` in general), which is what the protocol relies on.

## Features

- **`Quaternion`** — a quaternion ring element over `Z/pZ` with:
  - Addition, subtraction, and Hamilton multiplication (mod `p`)
  - Algebraic conjugate (`a − bi − cj − dk`) and multiplicative norm
  - Modular inverse via Fermat's little theorem (`p` must be prime)
  - Group-theoretic conjugation (`w⁻¹ · q · w`)
- **`PublicParams`** / **`Party`** — protocol setup and per-party state (private pattern + secret word, round-1/round-2 message logic)
- **`run_key_exchange`** / **`run_key_exchange_with_seeds`** — runs a complete two-party exchange and reports whether the derived keys matched
- **Session key derivation** — SHA3-256 over a fixed-width big-endian serialization of the shared quaternion
- **Utility helpers** — hex encoding/decoding, bit expansion, and Hamming distance between hex strings (handy for analyzing key randomness/divergence)

## Installation

Add the crate dependencies to your `Cargo.toml` (adjust the package name/path to wherever you vendor this crate):

```toml
[dependencies]
rand = "0.8"
rand_chacha = "0.3"
sha3 = "0.10"
```

Then either add this crate as a path/git dependency, or drop `aag_quaternion.rs` directly into your project as a module.

## Quick start

```rust
use rand_chacha::ChaCha8Rng;
use rand::SeedableRng;
use aag_quaternion::{PublicParams, run_key_exchange_with_seeds};

fn main() {
    // 1. Generate public parameters: prime p, 4 generators each for Alice/Bob,
    //    secret words of length 6.
    let mut rng_params = ChaCha8Rng::seed_from_u64(123);
    let params = PublicParams::generate(1019, 4, 4, 6, &mut rng_params);

    // 2. Run the exchange (seeded here for reproducibility; pass `None` for
    //    real entropy).
    let result = run_key_exchange_with_seeds(&params, Some(1), Some(2));

    assert!(result.keys_match);
    println!("Shared session key: {}", result.session_key_hex.unwrap());
}
```

For finer control over each party's randomness source (e.g. using `OsRng` or a CSPRNG of your choice instead of seeded `ChaCha8Rng`), use [`run_key_exchange`] directly:

```rust
use rand::rngs::OsRng;
use aag_quaternion::{PublicParams, run_key_exchange};

let mut rng_params = OsRng;
let params = PublicParams::generate(2_147_483_647, 6, 6, 8, &mut rng_params);

let mut rng_alice = OsRng;
let mut rng_bob = OsRng;
let result = run_key_exchange(&params, &mut rng_alice, &mut rng_bob);

assert!(result.keys_match);
```

## API overview

| Item | Purpose |
|---|---|
| `Quaternion::new(a, b, c, d, p)` | Construct a quaternion, reducing components mod `p` |
| `Quaternion::norm()` | Multiplicative norm `a² + b² + c² + d²` mod `p` |
| `Quaternion::inverse()` | Modular multiplicative inverse (requires prime `p`) |
| `Quaternion::conjugate_by(w)` | Group conjugation `w⁻¹ · q · w` |
| `random_invertible_quaternion(p, rng)` | Sample a random invertible quaternion |
| `product(&[Quaternion])` | Left-to-right product of a quaternion sequence |
| `PublicParams::generate(p, k, m, L, rng)` | Generate public generator sets and parameters |
| `Party::new(name, generators, L, rng)` | Create a party with a fresh private pattern/secret word |
| `Party::conjugate_set(...)` / `Party::derive_shared(...)` | Round-1 and round-2 protocol steps |
| `run_key_exchange(...)` / `run_key_exchange_with_seeds(...)` | Run a full exchange end-to-end |
| `derive_session_key(&Quaternion)` | SHA3-256-based key derivation from a shared quaternion |
| `hamming_distance_hex(h1, h2)` | Bit-level distance between two hex-encoded digests |

Full documentation with examples is available via `cargo doc --open`, since every public item has doc comments.

## Running tests

```bash
cargo test
```

The test suite checks:
- Associativity and correctness of quaternion multiplication/inversion
- Multiplicativity of the norm
- That `product()` matches manual left-to-right multiplication
- That the key exchange succeeds and produces matching, well-formed (64 hex char) session keys across many trials
- That different secret patterns produce different session keys

## Security status

This implementation is intended for **learning, experimentation, and protocol prototyping** — not production use. In particular:

- **No formal security analysis.** AAG-style protocols are only as strong as the underlying conjugacy/decomposition search problem in the chosen platform group. The hardness of these problems in a quaternion ring over `Z/pZ` (as opposed to braid groups, where AAG was originally proposed) has not been established here.
- **Small/toy moduli in tests.** The test suite uses `p = 1019` purely for fast, deterministic testing — this is far too small for any real security margin.
- **No side-channel hardening.** Arithmetic is not constant-time.
- **No authentication.** Like baseline Diffie–Hellman, this exchange alone doesn't authenticate either party and is vulnerable to man-in-the-middle attacks without an additional authentication layer.

If you want production-grade key exchange, use a well-reviewed, standardized primitive (e.g. X25519/X448, or a NIST-selected post-quantum KEM) via an established crate. Use this library to learn how non-commutative key exchange protocols are structured, not as a drop-in secure channel.

## License

Choose a license for your project (e.g. MIT or Apache-2.0) and add a `LICENSE` file — none is specified here by default.
