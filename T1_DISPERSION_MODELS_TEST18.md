# TEST 18 — T1-dispersion model fitting

Adds a model-fitting tool to the Relaxometry panel and to the T1 Dispersion Results window.

## Fitting domain
All dispersion models are fitted to R1 (s^-1) versus proton Larmor frequency. T1 display is obtained from 1000/R1 after fitting. This is intentional because independent relaxation contributions add in R1.

## Background models
1. **Power law (biological)**: `R1 = Rinf + A f^-beta`.
   This is a compact empirical tissue/macromolecular model and is often a useful description of broad biological NMRD profiles.
2. **BPP single tau (standard)**: one dipolar correlation time using the standard `J(w)+4J(2w)` spectral-density form.
3. **BPP two tau (biological)**: two dipolar correlation-time components to represent multiscale motion.
4. **Power law + BPP (biological)**: broad tissue power-law background plus one resolved dipolar dispersion.

## Optional 14N quadrupolar enhanced relaxation
The tool can add empirical positive Lorentzians in R1:
- high-field pair at 2.08 and 2.76 MHz (about 49 and 65 mT),
- all three tissue peaks at about 0.70, 2.08 and 2.76 MHz,
- or a single broad ~2 MHz QRE term.

The default Q Lorentzian HWHM is 0.20 MHz. The resonance centres are fixed because sparse imaging dispersion curves generally cannot robustly determine centre, width and amplitude simultaneously. A shared width can be fitted when enough field samples are available.

The reported Q-area for each Lorentzian is `pi * amplitude * HWHM`; the total Q-area is the sum.

## Outputs
The fitting window shows:
- log-log measured and fitted T1 or R1 dispersion,
- R1 residuals,
- fitted parameters,
- R^2, RMSE, AICc, BIC and residual degrees of freedom.

The most recent model fit is stored at:
- `scan.T1FitResults.dispersionModelFit`
- `scan.param.T1DispersionModelFit`

CSV/XLSX/workspace export is available from the fitting window.

## References used for model choices
- Classical BPP/dipolar spectral density: `J(w)=tau/(1+w^2 tau^2)`.
- Biological tissue NMRD profiles commonly exhibit power-law backgrounds.
- Tissue 14N quadrupolar-enhanced relaxation is observed near ~0.7, 2.1 and 2.8 MHz; the two high-field features correspond to approximately 49 and 65 mT.
- QRE can be represented empirically by Lorentzian terms added to the background R1 dispersion when the full microscopic quadrupolar model is not identifiable from sparse data.

The QRE implementation here is deliberately an empirical imaging-data model, not a full Fries/Belorizky microscopic quadrupolar model.
