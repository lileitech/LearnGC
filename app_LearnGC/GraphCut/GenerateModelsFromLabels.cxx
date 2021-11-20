//
// GenerateModelsFromLabels
//   Usage: GenerateModelsFromLabels InputVolume Startlabel Endlabel
//          where
//          InputVolume is a meta file containing a 3 volume of
//            discrete labels.
//          StartLabel is the first label to be processed
//          EndLabel is the last label to be processed
//          NOTE: There can be gaps in the labeling. If a label does
//          not exist in the volume, it will be skipped.
//
//
#include <vtkMetaImageReader.h>
#include <vtkImageAccumulate.h>
#include <vtkDiscreteMarchingCubes.h>
#include <vtkWindowedSincPolyDataFilter.h>
#include <vtkMaskFields.h>
#include <vtkThreshold.h>
#include <vtkGeometryFilter.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkSmartPointer.h>
#include <vtkSTLWriter.h>
#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPolyDataMapper.h>  
#include <vtkActor.h>  
#include <vtkRenderer.h>  
#include <vtkRenderWindow.h>  
#include <vtkRenderWindowInteractor.h>  
#include <sstream> // string to number conversion
#include<iostream>



int main (/*int argc, char *argv[]*/)
{
  // Create all of the classes we will need
	vtkSmartPointer<vtkMetaImageReader> reader =
		vtkSmartPointer<vtkMetaImageReader>::New();

  vtkSmartPointer<vtkImageAccumulate> histogram =
    vtkSmartPointer<vtkImageAccumulate>::New();
  vtkSmartPointer<vtkDiscreteMarchingCubes> discreteCubes =
    vtkSmartPointer<vtkDiscreteMarchingCubes>::New(); //marching cubes是经过处理之后的二值图，discrete marching cubes是没有经过处理的原始label


  vtkSmartPointer<vtkWindowedSincPolyDataFilter> smoother =
    vtkSmartPointer<vtkWindowedSincPolyDataFilter>::New();
  vtkSmartPointer<vtkThreshold> selector =
    vtkSmartPointer<vtkThreshold>::New();
  vtkSmartPointer<vtkMaskFields> scalarsOff =
    vtkSmartPointer<vtkMaskFields>::New();
  vtkSmartPointer<vtkGeometryFilter> geometry =
    vtkSmartPointer<vtkGeometryFilter>::New();
  vtkSmartPointer<vtkXMLPolyDataWriter> writer =
    vtkSmartPointer<vtkXMLPolyDataWriter>::New();

  // Define all of the variables
  unsigned int startLabel = 0;
  unsigned int endLabel = 72;
  unsigned int smoothingIterations = 15;
  double passBand = 0.001;
  double featureAngle = 120.0;

  // Generate models from labels
  // 1) Read the meta file
  // 2) Generate a histogram of the labels
  // 3) Generate models from the labeled volume
  // 4) Smooth the models
  // 5) Output each model into a separate file

  //读入label图像
  reader->SetFileName("E:/LA_Segmentation/TraningData/P01025_6/111.mha");
  /*reader->SetFileName("E:\\LA_Segmentation\\TraningData\\P01025_6\\LA_420.mha");*/
  reader->Update();

 /* histogram->SetInputConnection(reader->GetOutputPort());
  histogram->SetComponentExtent(0, endLabel, 0,0, 0,0);
  histogram->SetComponentOrigin(0, 0, 0);
  histogram->SetComponentSpacing(1, 1, 1);
  histogram->Update();*/

  //生成mesh
  discreteCubes->SetInputConnection(reader->GetOutputPort());
  discreteCubes->GenerateValues(
    endLabel - startLabel + 1, startLabel, endLabel);
 /* discreteCubes->SetSampleSpacing();*/
  discreteCubes->Update();

  //平滑
  smoother->SetInputConnection(discreteCubes->GetOutputPort());
  smoother->SetNumberOfIterations(smoothingIterations);
  smoother->BoundarySmoothingOff();
  smoother->FeatureEdgeSmoothingOff();
  smoother->SetFeatureAngle(featureAngle);
  smoother->SetPassBand(passBand);
  smoother->NonManifoldSmoothingOn();
  smoother->NormalizeCoordinatesOn();
  smoother->Update();

  //输出一些值，来检验mesh
  vtkIdList* list=vtkIdList::New();
  std::cout << smoother->GetOutput()->GetPoint(10)[0] << endl;//输入下标，得到下标点对应的坐标值
  std::cout << smoother->GetOutput()->GetPoint(10)[1] << endl;
  std::cout << smoother->GetOutput()->GetPoint(10)[2] << endl;

 smoother->GetOutput()->GetCellPoints(10,list);//输入下标，得到其周围（mesh）上的三角点
 for (size_t i = 0; i < list->GetNumberOfIds(); i++)
 {
	 std::cout << smoother->GetOutput()->GetPoint(list->GetId(i))[0] << endl;
	 std::cout << smoother->GetOutput()->GetPoint(list->GetId(i))[1] << endl;
	 std::cout << smoother->GetOutput()->GetPoint(list->GetId(i))[2] << endl;
 }

  selector->SetInputConnection(smoother->GetOutputPort());
  selector->SetInputArrayToProcess(0, 0, 0,
                                   vtkDataObject::FIELD_ASSOCIATION_CELLS,
                                   vtkDataSetAttributes::SCALARS);

  // Strip the scalars from the output
  scalarsOff->SetInputConnection(selector->GetOutputPort());
  scalarsOff->CopyAttributeOff(vtkMaskFields::POINT_DATA,
                               vtkDataSetAttributes::SCALARS);
  scalarsOff->CopyAttributeOff(vtkMaskFields::CELL_DATA,
                               vtkDataSetAttributes::SCALARS);

  geometry->SetInputConnection(scalarsOff->GetOutputPort());


  //for (unsigned int i = startLabel; i <= endLabel; i++)
  //{
  //  // see if the label exists, if not skip it
  //  double frequency =
  //    histogram->GetOutput()->GetPointData()->GetScalars()->GetTuple1(i);
  //  if (frequency == 0.0)
  //  {
  //    continue;
  //  }
  //  // select the cells for a given label
  //  selector->ThresholdBetween(i, i);


	vtkSmartPointer<vtkPolyDataMapper> mapper =
		vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper->SetInputData(smoother->GetOutput());

	vtkSmartPointer<vtkActor>  actor =
		vtkSmartPointer<vtkActor>::New();
	actor->SetMapper(mapper);

	vtkSmartPointer<vtkRenderer> render =
		vtkSmartPointer<vtkRenderer>::New();
	render->AddActor(actor);
	render->SetBackground(10, 0, 0);

	vtkSmartPointer<vtkRenderWindow> rw =
		vtkSmartPointer<vtkRenderWindow>::New();
	rw->AddRenderer(render);
	rw->SetSize(640, 480);
	rw->SetWindowName("PolyData Structure Learning");
	rw->Render();

	vtkSmartPointer<vtkRenderWindowInteractor> rwi =
		vtkSmartPointer<vtkRenderWindowInteractor>::New();
	rwi->SetRenderWindow(rw);
	rwi->Initialize();
	rwi->Start();
 //}
  return 0;
}
