# qmc_mt — Quantum Microtubule Coherence Framework

Repositorio de reproducción para el manuscrito de coherencia cuántica en microtúbulos.

## Estado de reproducción

- **Entrypoint canónico**: `reproduce_paper_results.py`
- **Wrapper de compatibilidad**: `reproduce_paper.py` (redirige al canónico)
- **Runbook release**: `RELEASE_RUNBOOK.md`
- **Artefactos principales**:
  - `validation_report.json` (auditable por máquina)
  - `LIVING_SI.md` (SI autogenerado)

## Ejecución rápida

Desde la carpeta `git_repo/`:

```bash
python reproduce_paper_results.py --fast
```

Ejecución completa por defecto:

```bash
python reproduce_paper_results.py
```

Barrido ROC extendido:

```bash
python reproduce_paper_results.py --full-roc
```

## Figuras del manuscrito

```bash
python generate_paper_figures.py
```

Las figuras canónicas se generan en `git_repo/figures_final/`.

## CI / calidad continua

El repo incluye workflow de CI en:

- `.github/workflows/ci.yml`

Cobertura mínima del workflow:

1. instalación editable del paquete,
2. tests unitarios (`tests/`),
3. test de consistencia NS (`test_ns_consistency.py`),
4. smoke Bayes HEOM v2 (`src/bayesian_heom_hierarchy_v2.py`),
5. smoke reproducible (`reproduce_paper_results.py --fast`),
6. smoke de figuras (`generate_paper_figures.py`).

## Jerarquía Bayesiana HEOM v2 (N pequeño)

- Script canónico: `src/bayesian_heom_hierarchy_v2.py`
- Input canónico: `src/heom_bayes_input_current.csv`
- Output canónico: `heom_bayes_out_v2/`

La v2 modela contracción en escala log para observables de salto positivos y aplica shrinkage jerárquico sobre `log(r)`. Está diseñada para **resumir evidencia de convergencia existente** en grupos pequeños (incluyendo grupos con 2 puntos), no para reemplazar nuevas corridas HEOM completas.

## Estructura relevante

- `src/qmc_mt/`: implementación física/estadística principal.
- `reproduce_paper_results.py`: pipeline E2E de validación.
- `reproduce_paper.py`: compatibilidad retroactiva del comando histórico.
- `generate_paper_figures.py`: generación de figuras finales.
- `paper.tex`, `paper.md`, `paper.bib`: manuscrito y bibliografía.

## Notas de release

- El comando recomendado para el paper es **solo uno**: `python reproduce_paper_results.py`.
- Si existen automatizaciones antiguas que llamen `reproduce_paper.py`, seguirán funcionando.
- Antes de release, ejecutar checklist completo de `RELEASE_RUNBOOK.md`.

## Estado editorial

- El manuscrito de software (`paper.md`) y bibliografía (`paper.bib`) se mantienen en este repo.
- Auditoría de consistencia editorial disponible en `EDITORIAL_AUDIT.md`.
- Nota: en `paper/paper.tex` hay claves bibliográficas aún no presentes en `paper.bib`; revisar `EDITORIAL_AUDIT.md` antes de cierre de versión.

## Licencia

MIT (`LICENSE`).
