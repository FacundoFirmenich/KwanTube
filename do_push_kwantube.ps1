# KwanTube — CI Fix Push
Set-Location "c:\Users\User\3D Objects\biofisicaquantiqaCLINE\KwanTube"

git add .github/workflows/ci.yml
git add tests/test_bayesian_heom_hierarchy_v2_smoke.py
git commit -m "fix: correct broken CI paths after script restructure and restrict CI to main branch only"
git push origin main

Write-Host "Done. CI: https://github.com/FacundoFirmenich/KwanTube/actions" -ForegroundColor Cyan
