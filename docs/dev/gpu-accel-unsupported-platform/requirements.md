# Requirements

## Goal
Hide GPU status fields from "CUDA Toolkit" through "Recommended" when the GPU Acceleration (Linux) section is shown on an unsupported platform.

## User story
- As a user on a non-Linux platform, I only see the Platform status and do not see CUDA-specific status fields that do not apply.

## Non-goals
- No backend or API changes.
- No changes to GPU install or refresh behavior.

## Acceptance criteria
- When GPU status reports the platform as unsupported, the UI does not render the CUDA Toolkit, CUDA Support, GPU4PySCF, cuTENSOR, or Recommended rows.
- When the platform is supported (Linux), the current status grid remains unchanged.
- Existing notices/errors and action buttons keep their current behavior.

## Open questions
- None.
