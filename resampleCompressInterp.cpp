#include <iostream>
#include <vector>

#include <vtkCellLocator.h>
#include <vtkDataSetWriter.h>
#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkResampleToImage.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkStructuredPointsReader.h>
#include <vtkStructuredPoints.h>

#include <mgard/compress.hpp>
#include <zstd.h>

//#define DEBUG_BUILD 1
#ifdef DEBUG_BUILD
#define DEBUG(x) std::cout << x << std::endl;
#else
#define DEBUG(x) \
  do             \
  {              \
  } while (0)
#endif

// input coordinates of 4 points
// input fieldvalue for 4 points

// input coordinates of dedicated point
// get the interpolated calue of dedicated point

template <std::size_t N, typename Real>
mgard::CompressedDataset<N, Real>
quantizeEncode(const mgard::TensorMeshHierarchy<N, Real> &hierarchy, Real *const u,
         const Real tolerance) {
  const std::size_t ndof = hierarchy.ndof();
  mgard::pb::Header header;
  mgard::populate_defaults(header);
  hierarchy.populate(header);
  mgard::decompose(hierarchy, header, u);
  {
    mgard::pb::ErrorControl &e = *header.mutable_error_control();
    e.set_mode(mgard::pb::ErrorControl::ABSOLUTE);
    e.set_norm(mgard::pb::ErrorControl::S_NORM);
    e.set_s(0);
    e.set_tolerance(tolerance);
  }
  mgard::MemoryBuffer<unsigned char> quantized = mgard::quantization_buffer(header, ndof);
  mgard::quantize(hierarchy, header, Real(0), tolerance, u, quantized.data.get());
  mgard::MemoryBuffer<unsigned char> buffer =
      mgard::compress(header, quantized.data.get(), quantized.size);
  return mgard::CompressedDataset<N, Real>(hierarchy, header, 0, tolerance,
                                    buffer.data.release(), buffer.size);
}

template <std::size_t N, typename Real>
Real *const
dequantizeDecode(const mgard::CompressedDataset<N, Real> &compressed) {
  const std::size_t ndof = compressed.hierarchy.ndof();
  Real *const dequantized = new Real[ndof];
  mgard::MemoryBuffer<unsigned char> quantized =
      mgard::quantization_buffer(compressed.header, ndof);
  mgard::decompress(compressed.header, const_cast<void *>(compressed.data()),
             compressed.size(), quantized.data.get(), quantized.size);
  mgard::dequantize(compressed, quantized.data.get(), dequantized);
  
  return dequantized;
}

double Interpolate2D(std::vector<std::vector<double>> inputCoor, std::vector<double> fieldValues, double dedicatedPoint[3])
{
  if (inputCoor.size() != fieldValues.size())
  {
    throw std::runtime_error("fields value should equal to the input coordinates");
  }

  // the point sequence is bottom left corner 0 (0,0), bottom right corner 1 (1,0)
  // upper left corner 2 (0,1) , upper right corner 3 (1,1)
  // https://www.omnicalculator.com/math/bilinear-interpolation
  double tolerance = 0.000001;
  double x1 = inputCoor[0][0] - tolerance;
  double y1 = inputCoor[0][1] - tolerance;

  double x2 = inputCoor[3][0] + tolerance;
  double y2 = inputCoor[3][1] + tolerance;

  // look at the bound, make sure the targeted point is in the bounds
  // there are tolerance issues

  // dedicatedPoint[0] should within x1 and x2
  // dedicatedPoint[1] should within y1 and y2
  double tx = dedicatedPoint[0];
  double ty = dedicatedPoint[1];
  if ((tx < x1) && (tx > x2))
  {
    DEBUG(tx << " , " << x1 << " , " << x2);
    throw std::runtime_error("tx is out of bounds");
  }

  if ((ty < y1) && (ty > y2))
  {
    DEBUG(ty << " , " << y1 << " , " << y2);
    throw std::runtime_error("ty is out of bounds");
  }

  double q11 = fieldValues[0];
  double q21 = fieldValues[1];
  double q12 = fieldValues[2];
  double q22 = fieldValues[3];

  DEBUG("x1 y1 x2 y2 " << x1 << "," << y1 << "," << x2 << "," << y2);
  DEBUG("q11 q21 q12 q22 " << q11 << "," << q21 << "," << q12 << "," << q22);

  double r1 = q11 * ((x2 - tx) / (x2 - x1)) + q21 * ((tx - x1) / (x2 - x1));
  double r2 = q12 * ((x2 - tx) / (x2 - x1)) + q22 * ((tx - x1) / (x2 - x1));

  double p = r1 * ((y2 - ty) / (y2 - y1)) + r2 * ((ty - y1) / (y2 - y1));
  DEBUG("interpolated p " << p);

  return p;
}


int main(int argc, char *argv[]) {
  // Parse command line arguments
  if (argc != 7) {
    std::cerr << "Usage: " << argv[0] << " datasetDirSuffix"  << " filedName" << " Smaple rate" ;
    std::cerr << " tolerance " << " split ratio (data : residual)";
    std::cerr << " residual compression strategy (0: lossless | 1: quantization + lossless) | 2: lossy"<< std::endl;
    return EXIT_FAILURE;
  }

  std::string datasetDirSuffix = argv[1];
  std::string fieldName = argv[2];
//  std::string unstructuredMeshRaw = datasetDirSuffix + ".vtk";
  std::string unstructuredWithVarFile = datasetDirSuffix + "WithVar.vtk";

  long unsigned int resampleNum=std::stoi(argv[3]);

  std::cout << "sample rate is " << resampleNum << std::endl;

  // load the vtk file
  vtkSmartPointer<vtkUnstructuredGridReader> reader =
      vtkSmartPointer<vtkUnstructuredGridReader>::New();
  reader->SetFileName(unstructuredWithVarFile.c_str());
  reader->Update();

  // get the specific polydata and check the results
  vtkUnstructuredGrid *unsGridData = reader->GetOutput();

  // unsGridData->Print(std::cout);
  // get field array direactly and put it into the data set
  vtkPointData *pointData = unsGridData->GetPointData();

  auto pointDataArray = pointData->GetScalars(fieldName.c_str());

  // pointDataArray->Print(std::cout);

  // set the pointDataArray into the mgard
  long unsigned int N = pointDataArray->GetNumberOfTuples();

  const mgard::TensorMeshHierarchy<1, double> hierarchy({N});

  double *const u = (double *)(pointDataArray->GetVoidPointer(0));
  double orignalDataSize = static_cast<double>(N * sizeof(*u));
  std::cout << "orignal element #: " << N << ", size = " << orignalDataSize/1024/1024 << " MB\n";

  // Now we set the compression parameters. First we select the norm in which to
  // control the compression error. We choose from the family of supported norms
  // by setting a smoothness parameter `s`. `s = 0` corresponds to the `L²`
  // norm.
  const double s = 0;
  // Next we set the absolute error tolerance `τ`. The approximating dataset `ũ`
  // generated by MGARD will satisfy `‖u - ũ‖_{L²} ≤ τ`.
  const double tolerance = std::stof(argv[4]);;
  const double ratio_t = std::stof(argv[5]); 
  if (ratio_t<=0.0) {
    std::cout << "tolerance of data cannot be zeros\n";
    exit(-1);
  }
  const double tol_data = ratio_t * tolerance;
  const double tol_resi = tolerance * (1.0-ratio_t); 
  
  int option = std::stoi(argv[6]);

  const mgard::CompressedDataset<1, double> compressed =
      mgard::compress(hierarchy, u, s, tolerance);
  std::cout << "compressed ok" << std::endl;
  const mgard::DecompressedDataset<1, double> decompressed = mgard::decompress(compressed);
  
  double min_v=1e9, max_v=0;
  for (size_t i=0; i<N; i++) {
    min_v = (u[i]<min_v) ? u[i] : min_v;
    max_v = (u[i]>max_v) ? u[i] : max_v;
  }
  std::cout << "value range: [" << min_v << ", " << max_v << "]\n";
  
  // `compressed` contains the compressed data buffer. We can query its size in
  // bytes with the `size` member function.
  
  double l2_err = 0.0, linf_err=0.0;
  size_t cnt_nzr=0;
  for (size_t i=0; i<N; i++) {
    double diff = std::abs(u[i]-decompressed.data()[i]);
    l2_err += diff*diff; 
    linf_err = (diff>linf_err) ? diff : linf_err;
    cnt_nzr = (u[i]>0) ? (cnt_nzr+1) : cnt_nzr;
  }
  l2_err = std::sqrt(l2_err / N);

  // try to create a sample based on
  // refer to this
  // https://gitlab.kitware.com/vtk/vtk/-/blob/master/Filters/Core/Testing/Cxx/TestResampleToImage.cxx
  vtkNew<vtkResampleToImage> resample;
  resample->SetUseInputBounds(true);
  // 2d case, the value at the last dim is 1
  resample->SetSamplingDimensions(resampleNum, resampleNum, 1);
  // resample->SetInputConnection(reader->GetOutputPort());
  resample->SetInputDataObject(unsGridData);
  resample->Update();

  vtkImageData *resampledImage = resample->GetOutput();
  // resampledImage->Print(std::cout);

  double bounds[6];
  resampledImage->GetBounds(bounds);

  std::cout << "resampled bounds: " ;

  //for (int i=0; i<6; i++){
  //  std::cout << bounds[i] << ",";
  //}

  std::cout << std::endl;

  double xspan = (bounds[1] - bounds[0])/49.0;
  //check the coordinates to make sure if it is correct
  //TODO, store this coordinates and put them into the vector
  //also store the y dim for this
  //for(int i=0;i<50;i++){
   // std::cout << "index " << i << " xcoord " << bounds[0]+i*xspan << std::endl;
  //}

  // check the coordinates of each point to see how it is connected with the grid
  
  for(int i=0; i < resampledImage->GetNumberOfPoints();i++){
    double tempCoor[3];
    resampledImage->GetPoint(i,tempCoor);
    //if(i<100){
    //    std::cout << "id " << i << " coor " << tempCoor[0] << "," << tempCoor[1] << "," << tempCoor[2] << std::endl;
    //}
  }

  // output resampled data for double checking
/*
  vtkSmartPointer<vtkDataSetWriter> writer =
      vtkSmartPointer<vtkDataSetWriter>::New();
  std::string outputFileName = fileSuffix + std::string("Resample.vtk");

  writer->SetFileName(outputFileName.c_str());

  // get the specific polydata and check the results
  writer->SetInputData(resampledImage);
  // Optional - set the mode. The default is binary.
  // writer->SetDataModeToBinary();
  // writer->SetDataModeToAscii();
  writer->Write();
*/
  // put into the mgard to check

  auto reamplePointDataArray =
      resampledImage->GetPointData()->GetScalars(fieldName.c_str());

  // pointDataArray->Print(std::cout);

  // set the pointDataArray into the mgard
  long unsigned int rN = reamplePointDataArray->GetNumberOfTuples();
  std::cout << "reamplePointDataArray points number " << rN  << ", ~ " << (float)rN / (float)N << "X of original size" << std::endl;
  const mgard::TensorMeshHierarchy<2, double> rhierarchy({resampleNum,resampleNum});

  double *const ru = (double *)(reamplePointDataArray->GetVoidPointer(0));

  const mgard::CompressedDataset<2, double> rcompressed =
      mgard::compress(rhierarchy, ru, s, tol_data);
  std::cout << "compressed sampled data ok" << std::endl;
  const mgard::DecompressedDataset<2, double> rdecompressed = mgard::decompress(rcompressed);

  // `compressed` contains the compressed data buffer. We can query its size in
  // bytes with the `size` member function.
  std::cout << "compression ratio for sampled 2d data set is: "
            << orignalDataSize / rcompressed.size()
            << ", compressed size: " << rcompressed.size()
            << std::endl;

  double l2_r2d = 0.0, linf_r2d = 0.0;
  for (size_t i=0; i<rN; i++) {
    double diff = std::abs(ru[i]-rdecompressed.data()[i]);
    l2_r2d += diff*diff; 
    linf_r2d = (diff > linf_r2d) ? diff : linf_r2d; 
  }
  l2_r2d = std::sqrt(l2_r2d / rN);
  std::cout << "resampled 2D's l2_err: " << l2_r2d << ", linf_err: " << linf_r2d << "\n"; 
  // TODO, use the compression based on 2d image

  const mgard::TensorMeshHierarchy<1, double> rhierarchy1d({resampleNum*resampleNum});

  const mgard::CompressedDataset<1, double> rcompressed1d =
      mgard::compress(rhierarchy1d, ru, s, tol_data);
  std::cout << "compressed sampled data ok" << std::endl;

  // `compressed` contains the compressed data buffer. We can query its size in
  // bytes with the `size` member function.
  std::cout << "compression ratio for sampled data set 1d is: "
            << orignalDataSize / rcompressed1d.size()
            << ", compressed size: " << rcompressed1d.size()
            << std::endl;

  // 

  // build the cell locator based on stuctured points
  vtkNew<vtkCellLocator> cellLocator;
  cellLocator->SetDataSet(resampledImage);
  cellLocator->BuildLocator();

  double pointcoord[3];
  auto fieldArray = resampledImage->GetPointData()->GetArray(fieldName.c_str());
  DEBUG("check field");
  DEBUG(fieldArray->GetNumberOfTuples());
  DEBUG(fieldArray->GetNumberOfComponents());

  std::vector<double> residual(N, 0);

  size_t cnt_resi = 0;
  double *d2d_data = (double*)rdecompressed.data();
  std::vector<double> interpoMGR(N); 
  for (vtkIdType id = 0; id < N; id++)
  {

    unsGridData->GetPoint(id, pointcoord);

    // use this pointcoord to find the specific cell in the structured grid
    vtkIdType interpCellId = cellLocator->FindCell(pointcoord);
    vtkCell *interpCell = resampledImage->GetCell(interpCellId);

    int interpPointsNum = interpCell->GetNumberOfPoints();

    if (interpCellId == -1)
    {
      // do not find the cell id
      // use the closetPoint to find
      double closestPoint[3];
      vtkIdType closestCellId;
      vtkGenericCell *closestCell;
      int subid = 0;
      double dist = 0;
      cellLocator->FindClosestPoint(pointcoord, closestPoint, closestCellId, subid, dist);
      interpCellId = closestCellId;
      interpCell = resampledImage->GetCell(interpCellId);
      DEBUG("original is -1 new id is " << closestCellId);
    }

    DEBUG("id " << id << " x " << pointcoord[0] << " y " << pointcoord[1]
                << " z " << pointcoord[2] << " interpCellId is " << interpCellId
                << " interpPointsNum is " << interpPointsNum);

    // get the cell from the cellid from resampledImage and do the interpolation
    vtkPoints *pointsArray = interpCell->GetPoints();
    int pointNum = pointsArray->GetNumberOfPoints();
    DEBUG(pointNum);
    vtkIdList *pointsIds = interpCell->GetPointIds();

    // std::cout << *(pointsArray->GetPoint(i)) << ",";
    std::vector<std::vector<double>> inputCoor;
    inputCoor.clear();
    std::vector<double> fieldValues;
    fieldValues.clear();
    for (int i = 0; i < pointNum; i++)
    {
      double strucCoors[3];
      pointsArray->GetPoint(i, strucCoors);

      DEBUG("point id " << pointsIds->GetId(i) << " coord " << strucCoors[0] << " " << strucCoors[1] << " " << strucCoors[2]);
      //  what are the value of assocaited points?
      //  get the field value for specific points
      double *fieldValue = &d2d_data[pointsIds->GetId(i)];//reamplePointDataArray->GetTuple(pointsIds->GetId(i));
      DEBUG("check data " << *fieldValue);
      //  use the strucCoors to interpolate pointcoord
      inputCoor.push_back({strucCoors[0], strucCoors[1], strucCoors[2]});
      fieldValues.push_back(*fieldValue);
    }

    interpoMGR.at(id) = Interpolate2D(inputCoor, fieldValues, pointcoord);
    double diff = interpoMGR.at(id) - *pointDataArray->GetTuple(id); 
//    std::cout << "id " << id << " diff: " << diff<< "\n"; 
    if (std::abs(diff) > tol_data) { 
        residual.at(id) = diff - tol_data; 
        cnt_resi ++;
    }
  }   
  
  std::cout << "number of residuals to be saved: " << cnt_resi << " (" << ((double)cnt_resi) / ((double)N)*100.0 << "%)\n";
  double *residualRCT = (double *)malloc(sizeof(double)*N);
  size_t resi_cSize;
  // compress the residual
  switch (option) {
      case 1:
        { 
          const size_t cBuffSize = ZSTD_compressBound(N);
          unsigned char *const zstd_resi = new unsigned char[cBuffSize];
          resi_cSize = ZSTD_compress(zstd_resi, cBuffSize, (void *)residual.data(), N, 1);
          ZSTD_decompress(residualRCT, N * sizeof(double), zstd_resi, resi_cSize);
          break;
        }
      case 2:
        {
          const mgard::CompressedDataset<1, double> encoded_resi =
            quantizeEncode(hierarchy, residual.data(), tol_resi); 
          residualRCT = dequantizeDecode(encoded_resi);
          break;
        }
      case 3: 
        {
          const mgard::CompressedDataset<1, double> compressed_resi = 
            mgard::compress(hierarchy, residual.data(), s, tol_resi);
          const mgard::DecompressedDataset<1, double> decompressed_resi = mgard::decompress(compressed_resi);
          resi_cSize = compressed_resi.size();
          memcpy(residualRCT, decompressed_resi.data(), sizeof(double)*N);
          break;
        }
  }
  // verify the statisfication of error bound
  double l2_comb = 0.0, linf_comb = 0.0;
  for (size_t id=0; id<N; id++) {
    double combValue = interpoMGR.at(id) - residualRCT[id];
    double diff = std::abs(combValue - *pointDataArray->GetTuple(id));
    l2_comb += diff*diff;
    linf_comb = (diff > linf_comb) ? diff : linf_comb;
  }
  l2_comb = std::sqrt(l2_comb / N);
  std::cout << "residual compression ratio bytes: " << resi_cSize << ", compression ratio using tol (" << tol_resi <<"): " << (orignalDataSize) / ((double)resi_cSize) << "\n";
  std::cout << "combined compression ratio for 2D resampled data: " << (orignalDataSize) / ((double)resi_cSize + rcompressed.size()) << "\n";
  std::cout << "interpolated data's l2_err: " << l2_comb << ", linf_err: " << linf_comb << "\n";
  std::cout << "compression ratio for raw data set is: "
            << orignalDataSize / compressed.size()
            << ", compressed size: " << compressed.size()
            << std::endl;
  std::cout << "original 1D's l2_err: " << l2_err << ", linf_err: " << linf_err << "\n"; 

  return 0;
}
