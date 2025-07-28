#include <iostream>

// STEP I/O includes
#include <STEPControl_Reader.hxx>
#include <STEPCAFControl_Writer.hxx>
#include <STEPCAFControl_Reader.hxx>
#include <TopoDS_Shape.hxx>
#include <IFSelect_ReturnStatus.hxx>

// XDE includes for assemblies and attributes
#include <TDocStd_Document.hxx>
#include <XCAFApp_Application.hxx>
#include <XCAFDoc_DocumentTool.hxx>
#include <XCAFDoc_ShapeTool.hxx>
#include <XCAFDoc_ColorTool.hxx>
#include <XCAFDoc_LayerTool.hxx>
#include <XCAFDoc_DocumentTool.hxx>
#include <TDataStd_Name.hxx>
#include <TDataStd_Comment.hxx>
#include <TDataStd_NamedData.hxx>
#include <TDF_Label.hxx>
#include <TDF_LabelSequence.hxx>

// Transformation includes
#include <gp_Trsf.hxx>
#include <gp_Vec.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <TopoDS_Compound.hxx>
#include <BRep_Builder.hxx>

// Color includes
#include <Quantity_Color.hxx>
#include <Quantity_NameOfColor.hxx>

// Validation properties
#include <GProp_GProps.hxx>
#include <BRepGProp.hxx>

// Shape exploration
#include <TopExp_Explorer.hxx>
#include <TopAbs_ShapeEnum.hxx>

// Tutorial
#include <TDocStd_Application.hxx>
#include <BinXCAFDrivers.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <XCAFDoc_DocumentTool.hxx>
#include <XCAFDoc_ColorTool.hxx>
#include <XCAFDoc_ShapeTool.hxx>
#include <TDF_ChildIterator.hxx>
#include <TopoDS.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <TopExp.hxx>
#include <gp_Quaternion.hxx>

TopoDS_Shape BuildWheel(const double OD, const double W){
  return BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(-W/2, 0, 0), gp::DX()), OD/2, W);
}

TopoDS_Shape  BuildAxle(const double D, const double L){
  return BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(-L/2, 0, 0), gp::DX()), D/2, L);
}

TopoDS_Shape BuildWheelAxle(TopoDS_Shape& wheel, TopoDS_Shape& axle, const double L){
  TopoDS_Compound compound;
  BRep_Builder builder;
  builder.MakeCompound(compound);

  gp_Trsf wheelT_right;
  wheelT_right.SetTranslationPart(gp_Vec(L/2, 0, 0));

  gp_Trsf wheelT_left;
  gp_Quaternion qn(gp::DY(), M_PI);
  gp_Trsf R;
  R.SetRotation(qn);
  wheelT_left = wheelT_left.Inverted() * R;
  wheelT_left.SetTranslationPart(gp_Vec(-L/2, 0, 0));
  
  builder.Add(compound, wheel.Moved(wheelT_right));
  builder.Add(compound, wheel.Moved(wheelT_left));
  builder.Add(compound, axle);

  return compound;
}


TopoDS_Shape BuildChassis(const TopoDS_Shape& wheelAxle, const double CL){
  TopoDS_Compound compound;
  BRep_Builder builder;
  builder.MakeCompound(compound);
  
  gp_Trsf frontT, rearT;
  frontT.SetTranslation(gp_Vec(0, CL/2, 0));
  rearT.SetTranslation(gp_Vec(0, -CL/2, 0));
  
  builder.Add(compound, wheelAxle.Moved(frontT));
  builder.Add(compound, wheelAxle.Moved(rearT));
  return compound;
}

struct t_prototype{
  TopoDS_Shape shape;
  TDF_Label label;
};

struct t_wheelPrototype : public t_prototype{
  TopoDS_Face frontFace;
  TDF_Label frontFaceLabel;
};

bool WriteStep(Handle(TDocStd_Document) doc, const std::string& filename) {
  STEPCAFControl_Writer writer;
  writer.SetColorMode(Standard_True);
  writer.SetNameMode(Standard_True);
  writer.SetLayerMode(Standard_True);
  writer.SetPropsMode(Standard_True);
  
  if (writer.Transfer(doc, STEPControl_AsIs)) {
      return writer.Write(filename.c_str()) == IFSelect_RetDone;
  }
  return false;
}

void tutorial(){

  std::cout << "Tutorial Start" << std::endl;

  Handle(TDocStd_Application) app = new TDocStd_Application;
  BinXCAFDrivers::DefineFormat(app);
  Handle(TDocStd_Document) doc;
  app->NewDocument("BinXCAF", doc);

  Handle(XCAFDoc_ShapeTool)
    ST = XCAFDoc_DocumentTool::ShapeTool(doc->Main());

    
  Handle(XCAFDoc_ColorTool)
    CT = XCAFDoc_DocumentTool::ColorTool(doc->Main());


  const double OD = 400;
  const double W = 100;
  const double D = 50;
  const double L = 500;
  const double CL = 1000;

  t_wheelPrototype wheelProto;
  wheelProto.shape = BuildWheel(OD, W);
  wheelProto.label = ST->AddShape(wheelProto.shape, false);


  t_prototype axleProto;
  axleProto.shape = BuildAxle(D, L);
  axleProto.label = ST->AddShape(axleProto.shape, false);

  t_prototype wheelAxleProto;
  wheelAxleProto.shape = BuildWheelAxle(wheelProto.shape, axleProto.shape, L);
  wheelAxleProto.label = ST->AddShape(wheelAxleProto.shape, true);

  t_prototype chassisProto;
  chassisProto.shape = BuildChassis(wheelAxleProto.shape, CL);
  chassisProto.label = ST->AddShape(chassisProto.shape, true);

  TDataStd_Name::Set(wheelProto.label, "wheel");
  TDataStd_Name::Set(axleProto.label, "axle");
  TDataStd_Name::Set(wheelAxleProto.label, "wheelAxle");
  TDataStd_Name::Set(chassisProto.label, "chassis");

  // Add metadata to test our new metadata export functionality
  Handle(TDataStd_NamedData) wheelMetadata = TDataStd_NamedData::Set(wheelProto.label);
  wheelMetadata->SetString("COMPAS_metadata_1", "My Awesome metadata! 1");
  wheelMetadata->SetString("COMPAS_metadata_2", "My Awesome metadata! 2");

  Handle(TDataStd_NamedData) axleMetadata = TDataStd_NamedData::Set(axleProto.label);
  axleMetadata->SetString("COMPAS_axle_metadata_1", "My Awesome metadata! 1");
  axleMetadata->SetString("COMPAS_axle_metadata_2", "My Awesome metadata! 2");

  Handle(TDataStd_NamedData) wheelAxleMetadata = TDataStd_NamedData::Set(wheelAxleProto.label);
  wheelAxleMetadata->SetString("COMPAS_wheel_axle_metadata_1", "My Awesome metadata! 1");
  wheelAxleMetadata->SetString("COMPAS_wheel_axle_metadata_2", "My Awesome metadata! 2");

  Handle(TDataStd_NamedData) chassisMetadata = TDataStd_NamedData::Set(chassisProto.label);
  chassisMetadata->SetString("COMPAS_chassis_metadata_1", "My Awesome metadata! 1");
  chassisMetadata->SetString("COMPAS_chassis_metadata_2", "My Awesome metadata! 2");

  // Find actual component references (not just direct children)
  TDF_LabelSequence components;
  ST->GetComponents(chassisProto.label, components);
  for (int i = 1; i <= components.Length(); i++) {
      TDF_Label componentRef = components.Value(i);
      TDataStd_Name::Set(componentRef, "wheel_axle_instance");
  }
  // Set different colors on different components, not the same label twice
  CT->SetColor(wheelProto.label, Quantity_Color(1, 0, 0, Quantity_TOC_RGB), XCAFDoc_ColorGen);  // Red wheel
  CT->SetColor(axleProto.label, Quantity_Color(0, 1, 0, Quantity_TOC_RGB), XCAFDoc_ColorGen);   // Green axle
  CT->SetColor(wheelAxleProto.label, Quantity_Color(0, 0, 1, Quantity_TOC_RGB), XCAFDoc_ColorGen); // Blue assembly
  
  TopTools_IndexedMapOfShape allWheelFaces;
  TopExp::MapShapes(wheelProto.shape, TopAbs_FACE, allWheelFaces);

  wheelProto.frontFace = TopoDS::Face(allWheelFaces(2));
  wheelProto.frontFaceLabel = ST->AddSubShape(wheelProto.label, wheelProto.frontFace);
  CT->SetColor(wheelProto.frontFaceLabel, Quantity_Color(0, 0, 1, Quantity_TOC_RGB), XCAFDoc_ColorSurf);   // Green axle
 

  
  // PCDM_StoreStatus status = app->SaveAs(doc, "tutorial.xbf");
  // if (status != PCDM_SS_OK) {
  //   std::cout << "Failed to save document" << std::endl;
  //   return;
  // }

  if (!WriteStep(doc, "tutorial.step")) {
    std::cout << "Failed to write STEP file" << std::endl;
    return;
  }

  std::cout << "Tutorial End" << std::endl;
}


// Main function for standalone application
int main() {
    std::cout << "Running OCCT STEP Tutorial Example..." << std::endl;
    
    // Run the tutorial function
    tutorial();
    
    std::cout << "Tutorial completed successfully!" << std::endl;
    return 0;
}

