#ifndef MRVM_H
#define MRVM_H

#include <list>
#include "matrix.h"
#include <QString>
#include <QMap>
#include <QFile>
#include <QFileInfo>
#include <arrayfire.h>
#include <QXmlStreamReader>
#include "gdal_priv.h"
#include "cpl_conv.h" // for CPLMalloc()

#ifndef NDEBUG
#   define ASSERT(condition, message) \
  do { \
  if (! (condition)) { \
  std::cerr << "Assertion `" #condition "` failed in " << __FILE__ \
  << " line " << __LINE__ << ": " << message << std::endl; \
  std::exit(EXIT_FAILURE); \
  } \
  } while (false)
#else
#   define ASSERT(condition, message) do { } while (false)
#endif

class MLTOOLS_EXPORT MRVMItem
{
  
public:
  enum MRVMValueType
  {
    Real,
    Categorical
  };

  enum IOType
  {
    Input,
    Output
  };

public:
  MRVMItem(IOType type = IOType::Input, const QString& name = "");
  
  ~MRVMItem();

  QString name() const;

  void setName(const QString& name);
  
  void clearAllValues();

  const QList<QString>& trainingValuesAsString() const;

  virtual float* trainingValues(int row) = 0;

  virtual void setTrainingValuesAsString(const QList<QString>& trainingValues);

  virtual void setTrainingValues(int row, float*& values) = 0;

  const QList<QString>& forecastValuesAsString() const;

  virtual float* forecastValues(int row) = 0;

  virtual void setForecastValuesAsString(const QList<QString>& forecastValues);

  virtual void setForecastValues(int row, float*& values) = 0;

  const QList<QString>& forecastUncertaintyValuesAsString() const;

  virtual float* forecastUncertaintyValues(int row) = 0;

  virtual void setForecastUncertaintyValueAsString(const QList<QString>& forecastUncertaintyValuesAsString);

  virtual void setForecastUncertaintyValues(int row, float*& values) = 0;

  virtual void readXML(QXmlStreamReader & xmlReader);

  virtual void writeXML(QXmlStreamWriter& xmlWriter);

  virtual int columnCount() = 0;
  
  int numTrainingValues() const;

  int numForecastValues() const;

  const QMap<QString, QString>& properties() const;

  virtual void setProperties(const QMap<QString, QString>& properties);
  
  virtual MRVMValueType valueType() const = 0;
  
  virtual QString type() const = 0;

  virtual IOType ioType() const;

  virtual QString toString() const;

protected:
  IOType m_iotype;
  QString m_name;
  QMap<QString, QString> m_properties;
  QList<QString> m_trainingValuesAsString,  m_forecastValuesAsString, m_forecastUncertaintyValuesAsString;

};

class MLTOOLS_EXPORT RealMRVMItem : public MRVMItem
{
  
public:
  
  RealMRVMItem(IOType type = IOType::Input,const QString& name = "");
  
  ~RealMRVMItem();
  
  float* trainingValues(int row) override;


  virtual void setTrainingValuesAsString(const QList<QString>& trainingValues) override;

  virtual void setTrainingValues(int row, float*& values)  override;


  float* forecastValues(int row) override;

  virtual void setForecastValues(int index, float*& values) override;


  float* forecastUncertaintyValues(int row) override;

  virtual void setForecastUncertaintyValues(int index, float*& values) override;


  int columnCount() override;

  MRVMValueType valueType() const override;

  QString type() const override;


};

class MLTOOLS_EXPORT RealArrayMRVMItem : public MRVMItem
{

public:
  RealArrayMRVMItem(IOType type = IOType::Input, const QString& name="");

  ~RealArrayMRVMItem();

  virtual float* trainingValues(int row) override;

  virtual void setTrainingValuesAsString(const QList<QString>& trainingValues) override;

  virtual void setTrainingValues(int row, float*& values) override;

  virtual float* forecastValues(int row) override;

  virtual void setForecastValues(int row, float*& values) override;

  virtual float* forecastUncertaintyValues(int row) override;

  virtual void setForecastUncertaintyValues(int row, float*& values) override;

  virtual void readXML(QXmlStreamReader & xmlReader) override;

  virtual int columnCount() override;

  virtual MRVMValueType valueType() const override;

  virtual QString type() const override;

protected:
  int m_columnCount;

};

class MLTOOLS_EXPORT CategoricalMRVMItem : public MRVMItem
{

public:
  CategoricalMRVMItem(IOType type = IOType::Input, const QString& name = "");

  ~CategoricalMRVMItem();

  virtual float* trainingValues(int row) override;

  virtual void setTrainingValues(int row, float*& values) override;

  virtual float* forecastValues(int row) override;

  virtual void setForecastValues(int row, float*& values) override;

  virtual float* forecastUncertaintyValues(int row) override;

  virtual void setForecastUncertaintyValues(int row, float*& values) override;

  virtual void readXML(QXmlStreamReader & xmlReader) override;

  virtual void writeXML(QXmlStreamWriter& xmlWriter) override;

  virtual int columnCount() override;

  virtual MRVMValueType valueType() const override;

  virtual QString type() const override;

  QMap<QString, int> categories() const;

protected:
  QMap<QString,int> m_categories;
  QMap<int,int> m_categoriesByIndex;
  QMap<int,int> m_indexByCategory;

};

class MLTOOLS_EXPORT RealRaster : public RealArrayMRVMItem
{

public:

  RealRaster(IOType type = IOType::Input,const QString& name ="");

  ~RealRaster();

  virtual float* trainingValues(int row) override;

  virtual void setTrainingValuesAsString(const QList<QString>& trainingValues) override;

  virtual void setTrainingValues(int row, float*& values) override;

  virtual float* forecastValues(int row) override;

  virtual void setForecastValues(int row, float*& values) override;

  virtual void setForecastValuesAsString(const QList<QString>& forecastValues) override;

  virtual void setForecastUncertaintyValueAsString(const QList<QString>& forecastUncertaintyValuesAsString);

  virtual float* forecastUncertaintyValues(int row) override;

  virtual void setForecastUncertaintyValues(int row, float*& values) override;

  virtual void readXML(QXmlStreamReader & xmlReader) override;

  virtual QString type() const override;

private:
  void writeDataToRaster(const QString& filePath, float*& values);

  float* readDataFromRaster(const QString& filePath);

  void readRasterProperties();

  void createOutputRasters();

private:
  int m_xSize, m_ySize;
  int* m_validCell;
  float m_noData;
  GDALDriver* m_driver;
  double m_gcp[6];
  const char* m_wktproj;
};

//class MLTOOLS_EXPORT CategoricalRaster : public CategoricalMRVMItem
//{

//public:

//  CategoricalRaster(IOType type = IOType::Input,const QString& name = "");

//  ~CategoricalRaster();

//  virtual float* trainingValues(int row) override;

//  virtual void setTrainingValuesAsString(const QList<QString>& trainingValues) override;

//  virtual void setTrainingValues(int row, float*& values) override;

//  virtual float* forecastValues(int row) override;

//  virtual void setForecastValues(int row, float*& values) override;

//  virtual void setForecastValuesAsString(const QList<QString>& forecastValues) override;

//  virtual void setForecastUncertaintyValueAsString(const QList<QString>& forecastUncertaintyValuesAsString);

//  virtual float* forecastUncertaintyValues(int row) override;

//  virtual void setForecastUncertaintyValues(int row, float*& values) override;

//  virtual void readXML(QXmlStreamReader & xmlReader) override;

//  virtual QString type() const override;

//private:
//  void writeDataToRaster(const QString& filePath, float*& values);

//  float* readDataFromRaster(const QString& filePath);

//  void readRasterProperties();

//  void createOutputRasters();

//private:
//  int m_xSize, m_ySize;
//  int* m_validCell;
//  float m_noData;
//  int m_columnCount;
//  GDALDriver* m_driver;
//  double m_gcp[6];
//  QString m_wktproj;
//};


class MLTOOLS_EXPORT Kernel
{
public:
  enum KernelType
  {
    Gaussian,
    Laplace,
    Polynomial,
    HomogeneousPolynomail,
    Spline,
    Cauchy,
    Cubic,
    Distance,
    ThinPlateSpline,
    Bubble
  };

public:
  Kernel(KernelType type = Gaussian , double lengthScale = 1000.0 );

  ~Kernel();
  
  KernelType type() const;
  
  void setKernelType(KernelType type);
  
  double lengthScale() const;
  
  void setLengthScale(double lengthScale);
  
  double polynomialPower() const;
  
  void setPolynomialPower(double polynomialPower);

  bool useBias() const;

  void setUseBias(bool m_useBias);
  
  af::array calculateKernel(const af::array& x1, const af::array& x2);

  //private:
public:

  af::array calculateGaussianKernel(const af::array& x1, const af::array& x2);

  af::array calculateLaplaceKernel(const af::array& x1, const af::array& x2);

  af::array calculatePolynomialKernel(const af::array& x1, const af::array& x2);
  
  af::array calculateHomogeneousPolynomialKernel(const af::array& x1, const af::array& x2);

  af::array calculateSplineKernel(const af::array& x1, const af::array& x2);

  af::array calculateCauchyKernel(const af::array& x1, const af::array& x2);

  af::array calculateCubicKernel(const af::array& x1, const af::array& x2);

  af::array calculateDistanceKernel(const af::array& x1, const af::array& x2);

  af::array calculateThinPlateSplineKernel(const af::array& x1, const af::array& x2);

  af::array calculateBubbleKernel(const af::array& x1, const af::array& x2);

  af::array distanceSquared(const af::array& x , const af::array& y);

private:
  KernelType m_kernelType;
  double m_lengthScale, m_polynomailPower;
  bool m_useBias;
};

class MLTOOLS_EXPORT MRVM : public QObject
{

public:
  enum Mode
  {
    TrainingAndRegression,
    Training,
    Regression
  };
  
  Q_OBJECT;
  
public:
  MRVM(const QFileInfo& file);
  
  ~MRVM();
  
  QString name() const;
  
  int inputDimension() const;
  
  int targetDimension() const;
  
  int maxNumberOfIterations() const;
  
  void setMaxNumberOfIterations(int niters);

  bool verbose() const;

  void setVerbose(bool verbose);

  int numberOfIterations() const ;
  
  MRVM::Mode mode() const;

  bool converged () const;

  const Kernel& kernel() const;
  
  const QList<MRVMItem*>& inputItems() const;
  
  void addInputItem(MRVMItem* const inputItem);
  
  bool removeInputItem(MRVMItem* const inputItem);
  
  const QList<MRVMItem*>& outputItems() const;
  
  void addOutputItem(MRVMItem* const inputItem);
  
  bool removeOutputItem(MRVMItem* const inputItem);
  
  QString matrixOutputFile() const;
  
  void setMatrixOutputFile(const QString& matrixOutputFile);
  
  const af::array& usedRelevantVectors() const;

  const af::array& alpha() const;
  
  const af::array& invSigma() const;

  const af::array& omega() const;

  const af::array& mu() const;
  
  void saveProject();

  void start();

  void performTraining();

  void mrvm();

  void fmrvm();

  void performRegression();

  static bool gdalRegistered();

#if Q_DEBUG
#else
private:
#endif

  void readProject();

  void qaqc();
  
  MRVMItem* readMRVMItem(MRVMItem::IOType type,  QXmlStreamReader& reader);
  
  float* getInputMatrix(int& row, int& column, bool training = true);

  float* getTargetMatrix(int& row, int& column, bool training = true);

  void writeOutput();

  QPair<int,double> calculateDeltaL(const af::array& mask, const af::array& alpha,
                                    const af::array& alphaNew, const af::array& sPrime,
                                    const af::array& qPrime, const af::array& s, const af::array& q);



  af::array corrcov(const af::array& cov);
  
private:
  QList<MRVMItem*> m_inputItems, m_outputItems;
  af::array m_inputMatrix, m_targetMatrix , m_used , m_alpha , m_invSigma, m_omega , m_inputForecastMatrix , m_Mu , m_outputForecastMatrix , m_stdPrediction;
  int m_maxNumberOfIterations = 1000 , m_numberOfIterations ;
  float m_minChangeAlpha, m_maxChangeAlpha;
  QString m_matrixOutputFile;
  int N, V;
  QFileInfo m_file;
  QString m_name;
  MRVM::Mode m_mode;
  Kernel m_kernel;
  static bool s_gdalRegistered;
  bool m_converged;
  double m_tolerance = 0.01;
  int algmode = 0;
  bool m_verbose;
};



#endif // MRVM_H
