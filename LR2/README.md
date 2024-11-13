Start


Чтобы все работало нужно в task.json в: "-I", "path\\eigen-3.4.0" - поменять на свой путь 

и в  c_cpp_properties.json : "includePath": ["${workspaceFolder}/**","path\\eigen-3.4.0"] - на свой путь

```
g++ -O3 -mavx2 -ftree-vectorize -march=native -I "path-library/eigen-3.4.0" -o lr2.exe lr2.cpp
./lr2
```