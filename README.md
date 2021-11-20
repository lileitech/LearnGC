# LearnGC

## Overview
The repository contains the core codes of "[Atrial scar quantification via multi-scale CNN in the graph-cuts framework](https://www.sciencedirect.com/science/article/pii/S1361841519301355)".
The resposutory includes three folds:
### app_MedAI2018 fold
This fold includes the original C++ scripts to generate the multi-scale patches and the pre-processing code for LearnGC.
### MS-CNN fold
This fold includes the python code to train and test the LearnGC architecture.
### Matlab_script fold
This fold includes some pre-processing scripts employed in AtrialJSQnet, and some of these scripts aimed to use the generated C++ tools mentioned in the C++ script fold.

Besides, GenerateModelsFromLabels.cxx was used to convert 3D label into 3D mesh.

## Releated work
You may also be interested in following papers:
1. [AtrialJSQnet: A New Framework for Joint Segmentation and Quantification of Left Atrium and Scars Incorporating Spatial and Shape Information](https://www.sciencedirect.com/science/article/pii/S1361841521003480)
2. [Medical Image Analysis on Left Atrial LGE MRI for Atrial Fibrillation Studies: A Review](https://arxiv.org/pdf/2106.09862.pdf)
3. [AtrialGeneral: Domain Generalization for Left Atrial Segmentation of Multi-center LGE MRIs](https://link.springer.com/chapter/10.1007/978-3-030-87231-1_54)


## Cite
If this code is useful for you, please kindly cite this work via:

@article{journal/MedIA/li2020,  
  title={Atrial scar quantification via multi-scale {CNN} in the graph-cuts framework},  
  author={Li, Lei and Wu, Fuping and Yang, Guang and Xu, Lingchao and Wong, Tom and Mohiaddin, Raad and Firmin, David and Keegan, Jennifer and Zhuang, Xiahai},  
  journal={Medical image analysis},  
  volume={60},  
  pages={101595},   
  year={2020},    
  publisher={Elsevier}  
}


If you have any questions, you are always welcome to contact with lilei.sky@sjtu.edu.cn.


