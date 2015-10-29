#ifndef MRVM_H
#define MRVM_H

#include <list>
#include "matrix.h"
#include <QString>
#include <QMap>
#include <QFile>
#include <QFileInfo>
#include "gdal_priv.h"
#include "cpl_conv.h"
#include <arrayfire.h>
#include <QXmlStreamReader>

using namespace af;

class MLTOOLS_EXPORT MRVMItem
{
  
public:
  enum MRVMValueType
  {
    Real,
    Categorical
  };
  
public:
  MRVMItem(const QString& name = "");
  
  ~MRVMItem();
  
  void setName(const QString& name);
  
  QString name() const;
  
  virtual void addValue(const QString& value);

  void clearValues();

  const QList<QString>& values() const;

  virtual double* values(int index) = 0;

  virtual void setValues(int index, double*& values) = 0;
  
  virtual int columnCount() = 0;
  
  virtual void setProperties(const QMap<QString, QString>& properties);
  
  const QMap<QString, QString>& properties() const;
  
  virtual MRVMValueType mRVMValueType() const = 0;
  
  virtual QString type() const = 0;

protected:
  QString m_name;
  QMap<QString, QString> m_properties;
  QList<QString> m_values;
};

class MLTOOLS_EXPORT RealMRVMItem : public MRVMItem
{
  
public:
  
  RealMRVMItem(const QString& name = "");
  
  ~RealMRVMItem();
  
  virtual double* values(int index);

  virtual void setValues(int index, double*& values);
  
  int columnCount() override;
  
  MRVMValueType mRVMValueType() const override;

  virtual QString type() const;

};

class MLTOOLS_EXPORT RealArrayMRVMItem : public MRVMItem
{
  
public:
  
  RealArrayMRVMItem(const QString& name = "");
  
  ~RealArrayMRVMItem();
  
  virtual double* values(int index);
  
  virtual void setValues(int index, double*& values);
  
  virtual int columnCount() override;
  
  MRVMValueType mRVMValueType() const override;
  
  virtual QString type() const;

protected:
  int m_columnCount;
  
};

class MLTOOLS_EXPORT CategoricalMRVMItem : public MRVMItem
{
  
public:
  CategoricalMRVMItem(const QString& name = "");
  
  ~CategoricalMRVMItem();
  
  virtual double* values(int index);
  
  virtual void setValues(int index, double*& value);
  
  virtual int columnCount() override;
  
  MRVMValueType mRVMValueType() const override ;
  
  virtual QString type() const;

  const QMap<int,QString>& categories() const;
  
  void setCategories(const QMap<int,QString>& categories);
  
protected:
  QMap<int,QString> m_categories;
  
};

class MLTOOLS_EXPORT RealRaster : public RealArrayMRVMItem
{
  
public:
  
  RealRaster(const QString& name);
  
  ~RealRaster();
  
  void addValue(const QString& value) override;
  
  double* values(int index) override;

  void setValues(int index, double*& values) override;
  
  int columnCount() override;
  
  virtual QString type() const;

private:
  void readRasterProperties();
  
private:
  int m_xSize, m_ySize;
  int* m_validCell;
  double m_noData;
  int m_columnCount;
  GDALDriver* m_driver;
  double m_gcp[6];
  QString m_wktproj;
};


class MLTOOLS_EXPORT CategoricalRaster : public CategoricalMRVMItem
{

public:

  CategoricalRaster(const QString& name);
  
  ~CategoricalRaster();

  void addValue(const QString& value) override;

  double* values(int index) override;

  void setValues(int index, double*& value) override;

  int columnCount() override;

  virtual QString type() const;

private:
  void readRasterProperties();

private:
  int m_xSize, m_ySize;
  int* m_validCell;
  double m_noData;
  int m_columnCount;
  GDALDriver* m_driver;
  double m_gcp[6];
  QString m_wktproj;
};


class MLTOOLS_EXPORT Kernel
{
public:
  enum KernelType
  {
    Gaussian,
    Laplace,
    Polynomial,
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
  
  array calculateKernel(const array& x1, const array& x2);

//private:
public:

  array calculateGaussianKernel(const array& x1, const array& x2);

  array calculateLaplaceKernel(const array& x1, const array& x2);

  array calculatePolynomialKernel(const array& x1, const array& x2);

  array calculateSplineKernel(const array& x1, const array& x2);

  array calculateCauchyKernel(const array& x1, const array& x2);

  array calculateCubicKernel(const array& x1, const array& x2);

  array calculateDistanceKernel(const array& x1, const array& x2);

  array calculateThinPlateSplineKernel(const array& x1, const array& x2);

  array calculateBubbleKernel(const array& x1, const array& x2);

  array distanceSquared(const array& x , const array& y);

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
  
  int trainingItemsCount() const;

  int maxNumberOfIterations() const;
  
  void setMaxNumberOfIterations(int niters);
  
  MRVM::Mode mode() const;
  
  const Kernel& kernel() const;
  
  const QList<MRVMItem*>& inputItems() const;
  
  void addInputItem(MRVMItem* const inputItem);
  
  bool removeInputItem(MRVMItem* const inputItem);
  
  const QList<MRVMItem*>& outputItems() const;
  
  void addOutputItem(MRVMItem* const inputItem);
  
  bool removeOutputItem(MRVMItem* const inputItem);
  
  const array& usedRelevantVectors() const;
  
  void save();

  void start();

  void performTraining();

  void mrvm1();

  void mrvm2();

  void fmrvm1();

  void performRegression();

#if Q_DEBUG
#else
private:
#endif

  void readFile();
  
  MRVMItem* readMRVMItem(QXmlStreamReader& xmlReader) ;
  
  void writeMRVMItem(MRVMItem* const item, QXmlStreamWriter& xmlReader);
  
  float* getInputMatrix(int& row, int& column);

  float* getTargetMatrix(int& row, int& column);

  void setTargetMatrixRegressionValues(float*& values);

  QPair<int,double> calculateDeltaL(const array& mask, const array& alpha, 
                                    const array& alphaNew, const array& sPrime, 
                                    const array& qPrime, const array& s, const array& q);

  array corrcov(const array& cov);
  
private:
  QList<MRVMItem*> m_inputItems, m_outputItems;
  array m_inputMatrix, m_targetMatrix , m_used , m_alpha , invSigma, m_omega;
  int m_maxNumberOfIterations = 1000 , m_averageNumberOfIterations , m_trainingItemsCount;
  QFileInfo m_file;
  QString m_name;
  MRVM::Mode m_mode;
  Kernel m_kernel;
  static bool gdalRegistered;
  double m_tolerance = 0.01;
  int algmode = 1;
};



#endif // MRVM_H
