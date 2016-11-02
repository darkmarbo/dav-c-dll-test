# dav-c-dll-test
dav c++ dll test-exe

1 EIGEN文件夹是程序源码中所使用的矩阵库，几乎所有数据类型都使用它，所以程序项目中必须包含它
   在VS2010中的配置方法：
    1 在项目属性
         C/C++
 	    常规
 	    附加包含目录，添加eigen3目录即可
    2 项目中程序包含
        #include<Eigen/Dense>