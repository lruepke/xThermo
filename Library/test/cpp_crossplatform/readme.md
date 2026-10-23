Test how to use c++ code library on ios platform.

* cpp目录下存放的是关于c++的代码和动态库生成的东西，生成完了之后可以安装到当前目录的install文件夹

```
cmake -DCMAKE_TOOLCHAIN_FILE=./ios.cmake -DIOS_PLATFORM=OS -Bframewrok.ios -GXcode ..
```

* 编译ios app
```
cmake -DCMAKE_TOOLCHAIN_FILE=./ios.cmake -DIOS_PLATFORM=OS -Bapp.ios -GXcode ..
```