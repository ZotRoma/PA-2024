Start


Чтобы все работало нужно в task.json в: "-I", "path\\eigen-3.4.0" - поменять на свой путь 

и в  c_cpp_properties.json : "includePath": ["${workspaceFolder}/**","path\\eigen-3.4.0"] - на свой путь

```
g++.exe -fdiagnostics-color=always -fopenmp -g "path\lr3.cpp" -I path\eigen-3.4.0 -o "path\lr3.exe"
./lr3
```