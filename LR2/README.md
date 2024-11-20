Start


Чтобы все работало нужно в task.json в: "-I", "path\\eigen-3.4.0" - поменять на свой путь 

и в  c_cpp_properties.json : "includePath": ["${workspaceFolder}/**","path\\eigen-3.4.0"] - на свой путь

```
g++ -O3 -mavx2 -ftree-vectorize -march=native -I "path-library/eigen-3.4.0" -o lr2.exe lr2.cpp
./lr2
```

Код не меняется, меняется только флаги компиляции

Флаг -O3 включает агрессивные оптимизации, включая автоматическую векторизацию.
(Флаги для SIMD) Эти флаги указывают компилятору использовать конкретные инструкции SIMD, если ваш процессор их поддерживает.
SSE (Streaming SIMD Extensions) — инструкции для базовой векторизации (обычно для старых процессоров): -msse -msse2 -msse3 -mssse3 -msse4
AVX (Advanced Vector Extensions) — более новые и эффективные инструкции, доступные на большинстве современных процессоров: -mavx -mavx2