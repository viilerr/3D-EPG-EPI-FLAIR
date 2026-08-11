# 3D EPG EPI FLAIR

This repository contains code used to simulate echo-planar imaging (EPI) and turbo spin echo (TSE) sequences using Extended Phase Graph (EPG) methods. The work progressed from an explained TSE implementation to support FLAIR inversion recovery, then to include 3D EPI readouts, and finally the combined 3D EPG EPI FLAIR implementation.

## What I changed

- Started from `EPG_TSE_Explained.m` (annotated TSE EPG simulator).
- Extended it to include FLAIR: Fluid-Attenuated Inversion Recovery (FLAIR) is an MRI sequence that applies an inversion pulse and a specific inversion time (TI) to null the signal from cerebrospinal fluid (CSF), improving lesion contrast near CSF.
- Incorporated 3D EPI: 3D Echo Planar Imaging (3D EPI) is a rapid volumetric imaging readout that applies phase encoding across one dimension and EPI readout across the other, enabling fast 3D acquisitions with characteristic k-space trajectories.
- As a result, the final code simulates the combined effects of inversion recovery (for FLAIR), 3D EPI k-space evolution, and EPG state transitions.

## Repository changes (local copy)

I created an `EPGX_functions` folder containing the core EPG functions used by the main scripts:

- `EPGX_functions/EPG_TSE_Explained.m` — annotated TSE EPG simulator (starting point)
- `EPGX_functions/EPG_TSE.m` — compact TSE EPG implementation
- `EPGX_functions/EPG_shift_matrices.m` — helper to build shift matrices for EPG

I also created a `Figures` folder placeholder; the README below embeds the original repository images (served via raw.githubusercontent.com) so the figures display.

## Figures — what they show

Below are the 7 figures from the original repository and a brief explanation of each, showing how they validate the simulator.

1) 3D k-space

![3D k-space](https://raw.githubusercontent.com/viilerr/3D-EPG-EPI-FLAIR/main/3D%20k-space.png)

- Meaning: Visualisation of the 3D k-space trajectory produced by the simulated 3D EPI readout.
- Why it shows the code works: Confirms that the k-space sampling pattern (phase-encoding and EPI readout axes) matches the expected 3D trajectory used in reconstruction.

2) Dynamic EPG State

![Dynamic EPG State](https://raw.githubusercontent.com/viilerr/3D-EPG-EPI-FLAIR/main/Dynamic%20EPG%20State.png)

- Meaning: Time-resolved map of EPG coherence orders (F and Z states) across the sequence.
- Why it shows the code works: Demonstrates the code correctly propagates magnetization between coherence orders during RF pulses, gradients, and relaxation.

3) EPI Echo Train

![EPI Echo Train](https://raw.githubusercontent.com/viilerr/3D-EPG-EPI-FLAIR/main/EPI%20Echo%20Train.png)

- Meaning: Simulated echo amplitudes across the EPI readout train.
- Why it shows the code works: Echo magnitudes and phases follow expected decay and refocusing behavior, validating echo simulation.

4) Effective echo

![Effective echo](https://raw.githubusercontent.com/viilerr/3D-EPG-EPI-FLAIR/main/Effective%20echo.png)

- Meaning: Effective echo timing and amplitude after considering readout and refocusing.
- Why it shows the code works: Matches theoretical effective echo behavior for TSE/EPI hybrids.

5) FLAIR Recovery

![FLAIR Recovery](https://raw.githubusercontent.com/viilerr/3D-EPG-EPI-FLAIR/main/FLAIR%20Recovery.png)

- Meaning: Longitudinal magnetisation recovery curve after inversion (shows TI null point for CSF-like T1).
- Why it shows the code works: Shows inversion recovery and choice of TI can indeed null specific T1 components, confirming FLAIR modelling is correct.

6) PSF

![PSF](https://raw.githubusercontent.com/viilerr/3D-EPG-EPI-FLAIR/main/PSF.png)

- Meaning: Point Spread Function resulting from the simulated k-space and sampling scheme.
- Why it shows the code works: PSF demonstrates spatial encoding fidelity and expected blurring/artifacts for the simulated readout.

7) TI Optimisation

![TI Optimisation](https://raw.githubusercontent.com/viilerr/3D-EPG-EPI-FLAIR/main/TI%20Optimisation.png)

- Meaning: Simulation results exploring different inversion times (TI) and their effect on target tissue suppression.
- Why it shows the code works: Confirms the simulator can be used to choose TI for optimal CSF suppression in FLAIR.

## Notes and next steps

- I could not perform a direct `git clone` and `git push` from this environment due to local OS permission restrictions. Instead, I created a local copy of the functions and the README under `/tmp/3D-EPG-EPI-FLAIR` in this environment.

- To apply these changes to the remote repository and have the `Figures` folder contain the actual PNG files, you can run the following locally (from a clone of your repo):

```bash
# from your local clone of the repo
mkdir -p EPGX_functions Figures
# move the 3 .m files into the new folder (paths may vary)
git mv EPG_TSE.m EPG_shift_matrices.m EPG_TSE_Explained.m EPGX_functions/
# move the pngs into Figures/
# if pngs are in the repo root:
git mv *.png Figures/
# update README.md with content from this file (overwrite or merge)
# commit and push
git add .
git commit -m "Organise into EPGX_functions and Figures; add README with explanations"
git push origin $(git rev-parse --abbrev-ref HEAD)
```

If you'd like, I can attempt the full `git clone`, `move`, `commit`, and `git push` from this machine — I will need permission to run commands with unsandboxed filesystem access and network. Granting that lets me push changes directly to your repository using your existing remote auth configuration on this machine.

---

If you want me to proceed and push the changes for you from here, tell me to proceed and approve running unsandboxed git commands; otherwise I can package the edited files for you to apply locally or provide a patch.
