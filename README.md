# MACROSCOPIC-MICROSCOPIC法

## 概要

このプロジェクトは、原子核物理学におけるFRLDM+Shell+BCSを解くためのプログラム群です．

## 参考文献

[FRDM2012](https://www.sciencedirect.com/science/article/pii/S0092640X1600005X)
[Axial Asymmetry](https://iopscience.iop.org/article/10.1088/0031-8949/8/1-2/003)





## 各moduleの依存関係


```mermaid
graph TD;
    constant_mod.f90-->nucleus_mod.f90;
    nucleus_mod.f90-->grid_mod.f90;
    grid_mod.f90-->frldm_mod.f90;
    CG_method_mod.f90-->frldm_mod.f90;
    grid_mod.f90-->micro_constant_mod.f90;
    micro_constant_mod.f90-->sp_solver_mod.f90
    CG_method_mod.f90-->sp_solver_mod.f90
    sp_solver_mod.f90-->shell_bcs_mod.f90
    frldm_mod.f90-->main.f90
    shell_bcs_mod.f90-->main.f90

```
                
## 運用

1. [プロジェクトのフォークを作成]
2. [変更を加える]
3. [プルリクエストを送信]



## 更新履歴

* 2026/4/21: [README.md を作成]
* 2026/5/14: [依存関係を明記]
