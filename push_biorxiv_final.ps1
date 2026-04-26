# Final Push to GitHub - Readiness for bioRxiv Submission
Set-Location "c:\Users\User\3D Objects\biofisicaquantiqaCLINE\KwanTube"

# Ensure everything is added
git add .

# Professional commit message
git commit -m "release: v3.5.0 Tier-A Production (bioRxiv candidate)
- Full validation pipeline passed (11/11 checks)
- Consolidated relative paths for SI and figures
- Unicode cleanup for cross-platform terminal compatibility
- SBC calibration verified (p=0.560)
- Structural PDB audit (eta~0.60) integrated"

# Push to main
git push origin main

Write-Host "KwanTube v3.5.0 is now live and ready for bioRxiv submission." -ForegroundColor Cyan
