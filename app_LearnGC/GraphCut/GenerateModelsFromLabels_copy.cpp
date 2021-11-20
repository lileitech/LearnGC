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
#include <vtkDiscreteMarchingCubes.h>
#include <vtkWindowedSincPolyDataFilter.h>
#include <vtkMaskFields.h>
#include <vtkThreshold.h>
#include <vtkGeometryFilter.h>
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
#include <iostream>

#include <vtkImageShiftScale.h>
#include <vtkNIFTIImageReader.h>
#include <vtkImageAccumulate.h>


int main (/*int argc, char *argv[]*/)
{

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
  vtkSmartPointer<vtkImageAccumulate> histogram =
	  vtkSmartPointer<vtkImageAccumulate>::New();
  vtkSmartPointer<vtkSTLWriter> stlwriter =
	  vtkSmartPointer<vtkSTLWriter>::New();


  // Define all of the variables
  unsigned int startLabel = 420;
  unsigned int endLabel = 500;
  unsigned int smoothingIterations = 15;
  double passBand = 0.001;
  double featureAngle = 120.0;

  std::string mainfold, casename,LAname, mesh_name;
  mainfold = "E:\\LA_Segmentation\\TraningData\\";
  casename = "P01025_6";
  LAname = mainfold + casename + "LA_label_1.nii.gz";
  mesh_name = mainfold + casename + "LA_mesh_1.nii.gz";

  vtkSmartPointer<vtkNIFTIImageReader> reader =
	  vtkSmartPointer<vtkNIFTIImageReader>::New();
  reader->SetFileName(LAname.c_str());
  reader->Update();

  /*
  //load dimensions using GetDataExtent
  int _extent = reader->GetDataExtent;
  //load spacing values
  float ConstPixelSpacing = reader->GetPixelSpacing;

  //输入图像原始的世界坐标信息
  vtkSmartPointer<vtkImageShiftScale> shiftScale =
	  vtkSmartPointer<vtkImageShiftScale> ::New();
  shiftScale->SetScale(reader->GetRescaleSlope());
  shiftScale->SetShift(reader->GetRescaleOffset());
  shiftScale->SetInputConnection(reader->GetOutputPort());
  shiftScale->Update();
  */

  
  /*
  reader->GetOutput()->GetDimensions(dims);
  std::cout << "原图像维度：" << dims[0] << " " << dims[1] << " " << dims[2] << std::endl;
  reader->GetOutput()->GetOrigin(origin);
  std::cout << "原图像原点：" << origin[0] << " " << origin[1] << " " << origin[2] << std::endl;
  reader->GetOutput()->GetSpacing(spacing);
  std::cout << "原图像间隔：" << spacing[0] << " " << spacing[1] << " " << spacing[2] << std::endl;

  vtkSmartPointer<vtkImageChangeInformation> changer =
	  vtkSmartPointer<vtkImageChangeInformation>::New();
  changer->SetInputData(reader->GetOutput());
//changer->SetOutputOrigin(178.1, -179.3, -77.78);
  changer->SetOutputSpacing(1, 1, 1);
//changer->SetCenterImage(0);
  changer->Update();

  changer->GetOutput()->GetDimensions(dims);
  std::cout << "修改后图像维度：" << dims[0] << " " << dims[1] << " " << dims[2] << std::endl;
  changer->GetOutput()->GetOrigin(origin);
  std::cout << "修改后图像中心：" << origin[0] << " " << origin[1] << " " << origin[2] << std::endl;
  changer->GetOutput()->GetSpacing(spacing);
  std::cout << "修改后图像间隔：" << spacing[0] << " " << spacing[1] << " " << spacing[2] << std::endl;
  */

  /*
  //clean all soft-tissue from the label image
  selector->SetInputConnection(reader->GetOutputPort());
  selector->ThresholdByLower(420);
  selector->ReplaceOutOn();
  selector->SetOutValue(0);  // set all values above 420 to 0
  selector->Update();
  */
  histogram->SetInputConnection(reader->GetOutputPort());
  histogram->SetComponentExtent(0, endLabel, 0, 0, 0, 0);
  histogram->SetComponentOrigin(0, 0, 0);
  histogram->SetComponentSpacing(1, 1, 1);
  histogram->Update();

  //生成mesh
  discreteCubes->SetInputConnection(reader->GetOutputPort());
  discreteCubes->GenerateValues(
    endLabel - startLabel + 1, startLabel, endLabel);
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


  /*
  //输出一些值，来检验mesh
  vtkIdList* list=vtkIdList::New();
  std::cout << smoother->GetOutput()->GetPoint(10)[0] <<" ";//输入下标，得到下标点对应的坐标值
  std::cout << smoother->GetOutput()->GetPoint(10)[1] << " ";
  std::cout << smoother->GetOutput()->GetPoint(10)[2] << endl;

 smoother->GetOutput()->GetCellPoints(10,list);//输入下标，得到其周围（mesh）上的三角点
 for (size_t i = 0; i < list->GetNumberOfIds(); i++)
 {
	 std::cout << smoother->GetOutput()->GetPoint(list->GetId(i))[0] <<" ";
	 std::cout << smoother->GetOutput()->GetPoint(list->GetId(i))[1] << " ";
	 std::cout << smoother->GetOutput()->GetPoint(list->GetId(i))[2] << endl;
 }
 */

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

  stlwriter->SetInputConnection(geometry->GetOutputPort());


  for (unsigned int i = startLabel; i <= endLabel; i++)
  {
	  // see if the label exists, if not skip it
	  double frequency =
		  histogram->GetOutput()->GetPointData()->GetScalars()->GetTuple1(i);
	  if (frequency == 0.0)
	  {
		  continue;
	  }

	  // select the cells for a given label
	  selector->ThresholdBetween(i, i);

	  // output the polydata
	  std::stringstream ss;
	  ss << mesh_name << i << ".vtp";
	  std::cout << " writing " << ss.str() << std::endl;

	  stlwriter->SetFileTypeToASCII();
	  stlwriter->SetFileName(ss.str().c_str());
	  stlwriter->Write();
  }

	/*vtkSmartPointer<vtkPolyDataMapper> mapper =
		vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper->SetInputData(smoother->GetOutput());

	vtkSmartPointer<vtkActor>  actor =
		vtkSmartPointer<vtkActor>::New();
	actor->SetMapper(mapper);

	vtkSmartPointer<vtkRenderer> render =
		vtkSmartPointer<vtkRenderer>::New();
	render->AddActor(actor);
	render->SetBackground(0, 0, 0);

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
	rwi->Start();*/

  return 0;
}
