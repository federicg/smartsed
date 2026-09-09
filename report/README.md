# report/

`cuda_port_benchmark.tex` — technical report on the GPU (CUDA) port of the
deterministic solver:

- milestone status and the kernel-port inventory (physical stage → CUDA wrapper);
- build/test infrastructure (`-DENABLE_CUDA=ON/OFF` from one source tree, the
  `-steps N` benchmark mode, `run/benchmark_sweep.sh`);
- correctness validation against the CPU reference via the `[FP]` field
  fingerprints (agreement ~1e-7 per step with `--use_fast_math`);
- an initial CPU-vs-GPU speed-up sweep over problem size (crossover near
  3–4 × 10⁴ active cells; ~2.2–2.3× on the GPU above it, on an NVIDIA A40);
- known issues and next steps.

## Build

```sh
cd report
latexmk -pdf cuda_port_benchmark.tex     # or: pdflatex cuda_port_benchmark.tex
```

Needs `pgfplots`, `siunitx`, `booktabs`, `fvextra`. LaTeX build products
(`*.aux`, `*.log`, `*.pdf`, …) are git-ignored.
