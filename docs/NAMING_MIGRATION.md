# Software naming migration — 4 October 2026

The pre-migration source is GitHub commit `e824281a2cce2012a0ccd6bce2585bf2d8b3e36f`. All 364 `.R`, `.py` and `.ps1` files were renamed. Directory `analysis/99_DIRECTION_UPGRADE/` became `analysis/direction_analysis/`. Names now use workflow-phase numbers, inventory numbers and descriptive functions; launchers share their target program prefix where a unique target is known.

`CODE_FILENAME_MAP.csv` records every original and new script path and the original SHA256. `RUN_ORDER.csv` lists the current software catalogue. `LEGACY_RUN_ORDER.csv` preserves the previous abbreviated overview. `START_HERE.md` identifies the tested portable entry points.

Scientific formulas, frozen eligibility rules, patient/sample identifiers, data paths and generated result names are unchanged. References to renamed software are updated; some hard-coded software locations are resolved relative to the distributed code. Historical raw-data workflows are still partly Windows/frozen-input bound and are not newly validated end to end.

Original `evidence/` files and the original `MANIFEST_SHA256.csv` preserve their historical names and hashes. They describe the prior execution, not the renamed distribution. The current manifest is `CODE_FILE_MANIFEST.csv`. New syntax and reproduction checks are in `evidence/NAMING_MIGRATION_20261004/`.

Inventory numbering is not a full executable dependency graph. Do not run all historical revisions, diagnostic trials and launchers sequentially. Missing historical sources and packaged assets are not recreated by this migration.

The unchanged proposed first Zenodo version is `1.0.0`; this is a draft version label, not a GitHub release or a published DOI. The original creator is Lu, Rongji and original software licensing remains MIT. No patient-level inputs are added to this code repository.
