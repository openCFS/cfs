# Characteristic (impedance-matched) flow ↔ acoustics coupling — derivation note

Spec for the openCFS-side boundary condition and the preCICE exchange used to couple
an OpenFOAM (compressible NS) domain to the openCFS acoustic wave equation across a
surface where only acoustic radiation is present.

**Target formulation: `ACOU_POTENTIAL`** (chosen because it avoids the node/element
pressure-result ambiguity and gives the particle velocity for free as a post-processing
result). The pressure formulation is kept as an appendix.

Status: DRAFT for review. No code changed yet (this file only).

---

## 0. Conventions

- `n` = **outward** unit normal of the domain the BC is written for.
- Perturbations primed: `p'`, `u'`; `u_n' = u' · n`.
- Reference (mean) acoustic medium `ρ₀`, `c₀`. **The same `ρ₀ c₀` reference must be used on
  both sides** so the interface impedance matches and the coupling is reflection-free. On the
  OpenFOAM side use this constant reference, not the local instantaneous `ρ c`.
- **openCFS potential convention (confirmed in code, AcousticPDE.cc:3190, 3202, 3232):**
  - particle velocity `v' = −∇ψ`  ⇒  `u_n' = v'·n = −∂ψ/∂n`
  - pressure `p' = ρ₀ ∂ψ/∂t`  (from linearized momentum `ρ₀ v'_t = −∇p'` with `v'=−∇ψ`)
- Acoustic characteristic invariants on the surface:
  - outgoing (+n):  `w_out = p' + ρ₀ c₀ u_n'`
  - incoming:       `w_in  = p' − ρ₀ c₀ u_n'`
  - one-way wave in +n: `p' = ρ₀ c₀ u_n'` ⇒ `w_out = 2·p_outgoing`.

Exchange rule: **each participant writes its own `w_out` and reads the neighbour's `w_out`
as its `w_in`** (consistent because `n_CFS = −n_OF`).

---

## 1. openCFS acoustic weak form (potential), as coded

Non-mech-coupled, time-domain potential formulation (`mechAcouFactor = 1`,
CalcMechAcouFac line 3736). Same wave operator as pressure, unknown `ψ`:

```
∫ (1/c₀²) ψ_tt q dV  +  ∫ ∇ψ · ∇q dV  −  ∮_Γ (∂ψ/∂n) q dS  =  0
```

- `coeffM = 1/c₀²`, `coeffK = 1` (stiffness `∫∇ψ·∇q`), confirmed lines ~441–476.
- ABC `coeffDamp = mechAcouFactor/c₀ = 1/c₀`, assembled as a `BBInt` boundary
  `IdentityOperator` into the `DAMPING` matrix (lines ~1848–1898).

The natural boundary datum is `∂ψ/∂n = −u_n'` — i.e. the normal velocity, which is exactly
why `normalVelocity` is the native BC for the potential formulation.

---

## 2. The impedance / characteristic boundary condition (potential)

Total = incident + scattered, `ψ = ψ_inc + ψ_sc`; absorbing condition applied to the
**scattered** part `∂ψ_sc/∂n = −(1/c₀) ψ_sc,t`. With `ψ_sc = ψ − ψ_inc`:

```
∂ψ/∂n + (1/c₀) ψ_t  =  ∂ψ_inc/∂n + (1/c₀) ψ_inc,t  ≡  g_ψ          (Robin BC)
```

Sanity of the absorbing half: for an outgoing wave `p' = ρ₀c₀ u_n'`, with `p'=ρ₀ψ_t` and
`u_n'=−∂ψ/∂n` ⇒ `∂ψ/∂n = −(1/c₀)ψ_t`. ✓ matches the existing ABC (`g_ψ = 0` ⇒ pure ABC).

Incident wave travels inward (−n): `∂ψ_inc/∂n = (1/c₀) ψ_inc,t`, so

```
g_ψ = (2/c₀) ψ_inc,t
```

Express via the exchanged invariant. `p'_inc = ρ₀ ψ_inc,t` ⇒ `ψ_inc,t = p'_inc/ρ₀`, and the
incident pressure is the neighbour's outgoing wave, `p'_inc = w_in/2`:

```
g_ψ = (2/c₀)·(p'_inc/ρ₀) = (2/(ρ₀c₀))·(w_in/2) = w_in/(ρ₀c₀)
```

**Key result — the forcing is ALGEBRAIC in `w_in` (no time derivative needed):**

```
g_ψ = w_in / (ρ₀ c₀),     with  w_in := w_out(neighbour).
```

### Weak-form split

```
−∮_Γ (∂ψ/∂n) q dS  =  ∮_Γ (1/c₀) ψ_t q dS   −   ∮_Γ g_ψ q dS
                      └──── LHS damping ────┘   └─── RHS load ───┘
```

- **LHS damping** `∮_Γ (1/c₀) ψ_t q dS` — **the existing `absorbingBCs` term**, reused verbatim.
- **RHS load** `f = ∮_Γ g_ψ q dS = ∮_Γ (w_in/(ρ₀c₀)) q dS` — the new incident-wave forcing.

> **Answer to "Robin, not a velocity pin?": yes.** The forcing is an inhomogeneous Neumann
> surface load (a `LinearForm`/`BUIntegrator`, RHS), assembled on the **same** surface as the
> `BBInt` damping term. Their **sum is the discrete Robin/impedance condition** above. It is
> **not** a Dirichlet `fixedValue` constraint (that would over-determine and reflect — the
> earlier failure mode). `w_in = 0` ⇒ load vanishes ⇒ pure absorbing BC (built-in check).

---

## 3. Mapping the RHS load onto `normalVelocity`

For non-mech `ACOU_POTENTIAL`, `CalcMechAcouFacWithCoef` passes the value through unchanged
(`scalFactor = 1`, `exValue = coef`, no density factor — line 3750), and the integrator
(lines 2216–2233) assembles

```
load = ∮_Γ v_fed · q dS
```

Match to `f = ∮_Γ (w_in/(ρ₀c₀)) q dS`:

```
v_fed = ± w_in / (ρ₀ c₀)
```

- Magnitude `w_in/(ρ₀c₀)`; it equals `2·(p'_inc/(ρ₀c₀))` = twice the physical incident normal
  velocity — the standard scattered-field doubling (consistent with reusing the absorbing term
  on the total field `ψ`).
- **Overall sign `±`:** depends on the LinearForm assembly sign for `normalVelocity` (whether
  the RHS contribution is `+∮v q` or `−∮v q`) combined with the `v=−∇ψ` orientation. Lock this
  empirically with the one-way-wave test (§6.2): a single global sign flip, trivial to verify.
- **No time derivative** of `w_in` is required — a real advantage of the potential path.

---

## 4. preCICE exchange

- **CFS writes** `w_out_CFS = p' + ρ₀ c₀ u_n'`, with
  - `p' = ρ₀ ψ_t` (physical pressure post-proc / `ACOU_POTENTIAL_DERIV_1`, line 3438),
  - `u_n' = ACOU_NORMAL_VELOCITY` (`= −∂ψ/∂n`, post-proc line 3231–3245).
  New scalar post-processing result `w_out` (or compute the combination in the adapter from the
  two existing results).
- **CFS reads** `w_in` → fed to `normalVelocity` as `v_fed = ± w_in/(ρ₀c₀)`.
- Add whitelist entries in `convertResultNamesToCFS` (PreciceAdapter.cc:234) for the new
  scalar(s) — exchange the **invariant `w` itself** (algebraic; no rate, unlike the existing
  `PressureTemporalDerivative` path).
- The read field reaches the BC through `CoefFunctionGridNodalDefaultPrecice` (existing wiring).

---

## 5. OpenFOAM side (mirror)

At the coupling patch, outward patch normal `n_OF`, reference `ρ₀ c₀`:

```
w_out_OF = (p − p₀) + ρ₀ c₀ ( (U − U₀) · n_OF )
```

- `p − p₀`: perturbation pressure (mean subtracted; adapter already distinguishes absolute vs
  perturbation pressure).
- `(U − U₀)·n_OF`: perturbation normal velocity (subtract mean flow if present).
- Computed in the FF adapter writer (or in the `codedMixed` patch BC, which needs `p,U` anyway).
- Patch BC is non-reflecting: reconstruct `p, U` from its own outgoing invariant (interior) +
  the incoming invariant read from openCFS — `codedMixed`/`waveTransmissive`-style, **not**
  `fixedValue`.
- No temporal derivative needed for the openCFS-bound data (`w` exchanged directly).

---

## 6. Sanity checks (build into the 1D testbed)

1. `w_in = 0` ⇒ RHS load vanishes ⇒ pure absorbing BC; must reproduce the existing ABC result.
2. Single one-way wave incident on the surface ⇒ fixes the global sign of `v_fed`; reflection
   coefficient ≈ 0.
3. Energy: `P = (w_out² − w_in²)/(4 ρ₀ c₀)` splits cleanly outgoing/incoming.
4. Forward pass: a wave generated in OpenFOAM crosses into acoustics; the far-field ABC
   back-reflection returns into OpenFOAM with correct amplitude/phase.

---

## 7. Code placement (per review feedback)

- **No new `DefineCharacteristicCouplingIntegrators()` function.** Add the BC inside the
  existing `AcousticPDE::DefineSurfaceIntegrators()`, in the `absorbingBCs` section (after the
  ABC loop / near the impedance TODO at ~line 1908), as a self-contained block that reuses the
  `coeffDamp` BBInt assembly and adds the `normalVelocity` LinearForm.
- New encapsulated XML element `characteristicCouplingBC` (type `DT_CharacteristicCouplingBC`)
  under `bcsAndLoads` in `CFS_PDEacoustic.xsd`.

---

## 8. Sign / convention quick-reference

| symbol | definition | source |
|---|---|---|
| `n` | outward normal | convention |
| `v' = −∇ψ`, `u_n' = −∂ψ/∂n` | velocity from potential | AcousticPDE.cc:3190/3232 |
| `p' = ρ₀ ψ_t` | pressure from potential | momentum + `v=−∇ψ` |
| `w_out` | `p' + ρ₀c₀ u_n'` | §0 |
| `w_in` | `p' − ρ₀c₀ u_n'` = neighbour's `w_out` | §0 |
| `coeffDamp` | `1/c₀` (ABC, reused) | AcousticPDE.cc:1848 |
| `g_ψ` (incident source) | `w_in/(ρ₀c₀)` — algebraic | §2 |
| `v_fed` (normalVelocity) | `± w_in/(ρ₀c₀)` | §3 |

---

## Appendix A. Pressure formulation (not chosen; kept for reference)

`ACOU_PRESSURE`, density-free form `∫(1/c₀²)p'_tt q + ∫∇p'·∇q − ∮∂p'/∂n q`, ABC `coeffDamp=1/c₀`.
There the incident forcing is `f = ∮ (1/c₀) ∂w_in/∂t q dS` (**requires the rate** `∂w_in/∂t`),
injected via `normalAcceleration` (the `normalVelocity` integrator throws for time-domain
pressure, lines 2261–2263). Matching its assembly `∮ ρ₀ a_fed q` gives
`a_fed = ∂w_in/∂t /(ρ₀c₀)`. Downsides vs. potential: needs the temporal derivative on the
exchange, and the pressure result's node/element ambiguity — hence the potential path was chosen.
