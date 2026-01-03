# Design

## Summary
Add a UI conditional in the GPU status grid to hide CUDA-related rows when `gpuStatus.is_linux` is false. Keep the Platform row and other UI elements (notices, errors, actions) unchanged.

## UI behavior
- Always render the Platform row.
- Render CUDA Toolkit, CUDA Support, GPU4PySCF, cuTENSOR, and Recommended rows only when `gpuStatus?.is_linux` is true.

## Implementation notes
- Update `src/web/pages/SettingsPage.tsx` to wrap the five CUDA-related grid items in a conditional block.
- Avoid touching the install/refresh buttons or status messages.

## Test plan
- Manual: simulate unsupported platform by setting `gpuStatus.is_linux = false` and confirm only the Platform row remains visible while the rest are hidden.
- Manual: simulate supported platform by setting `gpuStatus.is_linux = true` and confirm all rows render as before.
### Commands
- `npm run format:check`
- `npm run test:build`
- `npm run test:python-build`
- `/Users/goodapple/miniforge3/envs/pyscf-env/bin/python -m pytest src/python/tests -v`
