#ifndef MRVM_H
#define MRVM_H

#include "mrvm_global.h"
#include <list>
#include <QString>
#include <QMap>
#include <QFile>
#include <QPoint>
#include <QPointF>
#include <QFileInfo>
#include <QVariant>
#include <QPolygon>
#include <QSet>
#include <arrayfire.h>
#include <QXmlStreamReader>
#include "gdal_priv.h"
#include "cpl_conv.h" // for CPLMalloc()
#include <QRegularExpression>
#include <random>

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


# ifndef NUM_THREADS
#define NUM_THREADS 4
#endif

static std::default_random_engine gen;

class MRVM_EXPORT MRVMItem
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
    MRVMItem(IOType type = MRVMItem::Input, const QString& name = "");

    virtual ~MRVMItem();

    QString name() const;

    void setName(const QString& name);

    void clearAllValues();

    const QList<QString>& trainingValuesAsString() const;

    virtual af::array trainingValues(int valueIndex, int startRow, int length ) = 0;

    virtual void setTrainingValuesAsString(const QList<QString>& trainingValues);

    const QList<QString>& forecastValuesAsString() const;

    virtual af::array forecastValues(int valueIndex, int startRow, int length) = 0;

    virtual void setForecastValuesAsString(const QList<QString>& forecastValues);

    virtual void setForecastValues(int valueIndex, const af::array& values, const af::array& uncertainty) = 0;

    const QList<QString>& forecastUncertaintyValuesAsString() const;

    virtual void setForecastUncertaintyValueAsString(const QList<QString>& forecastUncertaintyValuesAsString);

    virtual void readXML(QXmlStreamReader & xmlReader);

    virtual void writeXML(QXmlStreamWriter& xmlWriter);

    virtual int columnCount() = 0;

    virtual int numTrainingValues() const;

    virtual int numForecastValues() const;

    virtual int numRowsPerTrainingValue() const;

    virtual int numRowsPerForecastValue() const;

    const QMap<QString, QVariant>& properties() const;

    virtual void setProperties(const QMap<QString, QVariant>& properties);

    virtual MRVMValueType valueType() const = 0;

    virtual QString type() const = 0;

    virtual IOType ioType() const;

    virtual QString toString() const;

protected:
    IOType m_iotype;
    QString m_name;
    QMap<QString, QVariant> m_properties;
    QList<QString> m_trainingValuesAsString,  m_forecastValuesAsString, m_forecastUncertaintyValuesAsString;
    int m_numTrainingValues, m_numForecastValues;

};

class MRVM_EXPORT RealMRVMItem : public MRVMItem
{

public:

    RealMRVMItem(IOType type = MRVMItem::Input,const QString& name = "");

    ~RealMRVMItem();

    af::array trainingValues(int valueIndex, int startRow = 0, int length = 1) override;

    virtual void setTrainingValuesAsString(const QList<QString>& trainingValues ) override;

    af::array forecastValues(int valueIndex, int startRow = 0, int length = 1) override;

    virtual void setForecastValues(int valueIndex, const af::array& values, const af::array& uncertainty) override;

    virtual void readXML(QXmlStreamReader & xmlReader) override;

    virtual void writeXML(QXmlStreamWriter& xmlWriter) override;

    int columnCount() override;

    virtual int numTrainingValues() const override;

    virtual int numForecastValues() const override;

    MRVMValueType valueType() const override;

    QString type() const override;

private:
    void readWriteTrainingValuesFiles(bool read = true);

    void readWriteForecastValuesFiles(bool read = true);

    void readWriteForecastUncertaintyValuesFiles(bool read = true);

    void expandListTo(QList<float>& list, int index);

    void expandListTo(QList<QString>& list, int index);

private:
    QList<float> m_trainingValues;
    QList<float> m_forecastValues;
    QList<float> m_forecastUncertaintyValues;

};

class MRVM_EXPORT RealArrayMRVMItem : public MRVMItem
{

public:
    RealArrayMRVMItem(IOType type = MRVMItem::Input, const QString& name="");

    ~RealArrayMRVMItem();

    virtual af::array trainingValues(int valueIndex, int startRow = 0, int length = 1) override;

    virtual void setTrainingValuesAsString(const QList<QString>& trainingValues) override;

    virtual af::array forecastValues(int valueIndex, int startRow = 0, int length = 1) override;

    virtual void setForecastValues(int valueIndex, const af::array& values, const af::array& uncertainty) override;

    virtual void readXML(QXmlStreamReader & xmlReader) override;

    virtual void writeXML(QXmlStreamWriter& xmlWriter) override;

    virtual int columnCount() override;

    virtual int numTrainingValues() const override;

    virtual int numForecastValues() const override;

    virtual MRVMValueType valueType() const override;

    virtual QString type() const override;

private:
    void readWriteTrainingValuesFiles(bool read = true);

    void readWriteForecastValuesFiles(bool read = true);

    void readWriteForecastUncertaintyValuesFiles(bool read = true);

    void expandListTo(QList<QList<float> >& list, int index);

    void expandListTo(QList<QString>& list, int index);

private:
    QList<QList<float> > m_trainingValues;
    QList<QList<float> > m_forecastValues;
    QList<QList<float> > m_forecastUncertaintyValues;
    static QRegularExpression s_regex;

protected:
    int m_columnCount;

};

class MRVM_EXPORT CategoricalMRVMItem : public MRVMItem
{

public:
    CategoricalMRVMItem(IOType type = MRVMItem::Input, const QString& name = "");

    ~CategoricalMRVMItem();

    virtual QMap<QString,int> getCategories() const;

    void setCategories(const QMap<QString,int>& categories);

    virtual af::array trainingValues(int valueIndex, int startRow = 0, int length = 1) override;

    virtual af::array forecastValues(int valueIndex, int startRow = 0, int length = 1) override;

    virtual void setForecastValues(int row, const af::array& values,
                                   const af::array& uncertainty) override;

    virtual void readXML(QXmlStreamReader & xmlReader) override;

    virtual void writeXML(QXmlStreamWriter& xmlWriter) override;

    virtual int columnCount() override;


    virtual int numTrainingValues() const override;

    virtual int numForecastValues() const override;


    virtual MRVMValueType valueType() const override;

    virtual QString type() const override;

    virtual QString toString() const override;

    QMap<QString, int> categories() const;

private:
    void readWriteTrainingValuesFiles(bool read = true);

    void readWriteForecastValuesFiles(bool read = true);

    void readWriteForecastUncertaintyValuesFiles(bool read = true);

    void expandListTo(QList<int>& list, int index);

    void expandListTo(QList<float>& list, int index);

    void expandListTo(QList<QString>& list, int index);

private:
    QList<int> m_trainingValues;
    QList<int> m_forecastValues;
    QList<float> m_forecastUncertaintyValues;

protected:
    QMap<QString,int> m_classbycategory;
    QMap<int,QString> m_categorybyclass;
    QMap<int,int> m_classbyindex;
    QMap<int,int> m_indexbyclass;\
    float maxCValue, minCValue;

};

class RasterBootstrapSampler;

class MRVM_EXPORT RasterItem
{

    friend class RasterBootstrapSampler;

public:
    RasterItem();

    virtual ~RasterItem();

    virtual QString getName() const = 0;

    virtual void resetProperties() = 0;

    bool contains(const QPointF& point, QPoint& pointIndex);

    bool isValid(const QPoint& index);

    bool isValid(int i , int j);

    QPointF getCoordinates(const QPoint& indexes)  const;

    QPointF getCoordinates(int x , int y)  const;

    QPoint getCoordinateIndexes(const QPointF& coordinates) const;

    void setBootstrapSamplingPoints(const QList<QPointF>& sampleLocations);

    bool includeLocation() const;

    virtual QPolygonF boundary() const;

protected:
    QPolygonF m_boundary;
    bool m_useRasterBootstrap;
    bool m_includeLocation;
    QList<QPoint> m_sampleLocations;
    int m_xSize, m_ySize;
    int* m_validCell;
    float m_noData;
    GDALDriver* m_driver;
    double m_gcp[6];
    const char* m_wktproj;
    RasterBootstrapSampler* m_rasterBootstrapSampler;
};

class MRVM_EXPORT RealRaster : public RealArrayMRVMItem, public RasterItem
{

public:

    RealRaster(IOType type = MRVMItem::Input,const QString& name ="");

    ~RealRaster();

    QString getName() const override;

    virtual af::array trainingValues(int valueIndex, int startRow, int length) override;

    virtual void setTrainingValuesAsString(const QList<QString>& trainingValues) override;

    virtual af::array forecastValues(int valueIndex, int startRow, int length) override;

    virtual void setForecastValues(int valueIndex, const af::array& values, const af::array& uncertainty) override;

    virtual void setForecastValuesAsString(const QList<QString>& forecastValues) override;

    virtual void setForecastUncertaintyValueAsString(const QList<QString>& forecastUncertaintyValuesAsString) override;

    virtual void readXML(QXmlStreamReader & xmlReader) override;

    virtual int numRowsPerTrainingValue() const override;
    virtual int numRowsPerForecastValue() const override;

    virtual void resetProperties() override;

    virtual QString type() const override;

private:

    void writeDataToRaster(const QString& filePath, const af::array& values);

    af::array readDataFromRaster(const QString& filePath, int startRow, int length);

    af::array readTrainingDataFromSampler(const QString& filePath , int startRow, int length);

    af::array readForecastDataFromSampler(const QString& filePath , int startRow, int length);

    void readRasterProperties();

    void createOutputRasters();

protected:
    int m_numRowsPerTrainingValue;
    int m_numRowsPerForecastValue;
    int m_numValidPixels;


};

class MRVM_EXPORT CategoricalRaster : public CategoricalMRVMItem , public RasterItem
{

public:

    CategoricalRaster(IOType type = MRVMItem::Input,const QString& name = "");

    ~CategoricalRaster();

    QString getName() const override;

    virtual af::array trainingValues(int valueIndex, int startRow , int length ) override;

    virtual void setTrainingValuesAsString(const QList<QString>& trainingValues) override;

    virtual af::array forecastValues(int valueIndex, int startRow, int length) override;

    virtual void setForecastValues(int valueIndex, const af::array& values, const af::array& uncertainty) override;

    virtual void setForecastValuesAsString(const QList<QString>& forecastValues) override;

    virtual void setForecastUncertaintyValueAsString(const QList<QString>& forecastUncertaintyValuesAsString);

    virtual void readXML(QXmlStreamReader & xmlReader) override;

    virtual int numRowsPerTrainingValue() const override;

    virtual int numRowsPerForecastValue() const override;

    virtual int columnCount() override;

    virtual void resetProperties() override;

    virtual QString type() const override;

private:

    void writeDataToRaster(const QString& filePathForecast, const af::array& values, const QString& filePathUncertain, const af::array& uncert);

    af::array readDataFromRaster(const QString& filePath, int startRow , int length);

    af::array readTrainingDataFromSampler(const QString& filePath, int startRow , int length);

    af::array readForecastDataFromSampler(const QString& filePath, int startRow , int length);

    void readRasterProperties();

    void createOutputRasters();

private:
    int m_numRowsPerTrainingValue;
    int m_numRowsPerForecastValue;
    int m_columnCount;
    int m_numValidPixels;

};

class MRVM_EXPORT RasterBootstrapSampler
{

public:
    RasterBootstrapSampler();

    ~RasterBootstrapSampler();

    int numSamples() const;

    void setNumSamples(int numSamples);

    QList<QString> rasterItemNames() const;

    QMap<QString, RasterItem*> rasterItems() const;

    void addRasterItem(RasterItem* rasterItem);

    bool removeRasterItem(RasterItem* rasteriterm);
    
    QList<QPointF> samplingLocations() const;

    void readXML(QXmlStreamReader & xmlReader);

    void writeXML(QXmlStreamWriter& xmlWriter);

    void createValidWindowCenters();

    void setRasterItemSamplingAttributes();

private:

    int m_numSamples;
    QList<QString> m_rasterItemNames;
    QMap<QString,RasterItem*> m_rasterItems;
    QList<QPointF> m_samplingLocations;

};

class MRVM_EXPORT Kernel
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

class MRVM_EXPORT MRVM : public QObject
{
public:
    enum Mode
    {
        TrainingAndRegression,
        Training,
        Regression
    };

public:
    MRVM(const QFileInfo& file);

    ~MRVM();

    QString name() const;

    int maxNumberOfIterations() const;

    void setMaxNumberOfIterations(int niters);

    bool verbose() const;

    void setVerbose(bool verbose);

    int numberOfIterations() const ;

    MRVM::Mode mode() const;

    bool converged() const;

    const Kernel& kernel() const;

    const QMap<QString, MRVMItem*>& inputItems() const;

    void addInputItem(MRVMItem*& inputItem);

    bool removeInputItem(const QString& name);

    const QMap<QString, MRVMItem*>& outputItems() const;

    void addOutputItem(MRVMItem*& inputItem);

    bool removeOutputItem(const QString& name);

    QList<RasterBootstrapSampler*> rasterBootstrapSamplers() const;

    QString matrixOutputFile() const;

    void setMatrixOutputFile(const QString& matrixOutputFile);

    const af::array& inputMatrix() const;

    const af::array& targetMatrix() const;

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

private:

    void readProject();

    void processRasterBootstrapSamplers();

    MRVMItem* readMRVMItem(MRVMItem::IOType type,  QXmlStreamReader& reader);

    void validateInputs();

    int ioValuesColumnCount(bool input = true) const;

    int ioValuesRowCount(int& maxrowsPerValue, bool input = true, bool training = true) const;

    af::array getInputMatrix(bool training = true);

    af::array getInputMatrix(int valueIndex, bool training = true);

    af::array getInputMatrix(int valueIndex, int startRow = 0, int length = 1, bool training = true);

    af::array getOutputMatrix(bool training = true);

    af::array getOutputMatrix(int valueIndex, bool training = true);

    af::array getOutputMatrix(int valueIndex,int startRow = 0, int length = 1, bool training = true);

    void writeOutput(int valueIndex, const af::array& values , const af::array& uncertainty);

    QPair<int,double> calculateDeltaL(const af::array& mask, const af::array& alpha,
                                      const af::array& alphaNew, const af::array& sPrime,
                                      const af::array& qPrime, const af::array& s,
                                      const af::array& q);

    af::array corrcov(const af::array& cov);

private:
    QMap<QString,MRVMItem*> m_inputItems, m_outputItems;
    QList<RasterBootstrapSampler*> m_rasterBootstrapSamplers;
    af::array m_inputMatrix, m_targetMatrix , m_used , m_alpha , m_invSigma, m_omega, m_Mu ;

    int N,V, m_maxNumberOfIterations = 1000 , m_numberOfIterations,
        m_numInputCols, m_numOutputCols, m_numInputTrainingRows, m_numOutputTrainingRows,
        m_numInputForecastRows, m_numOutputForecastRows, m_maxNumRowsPerForecastValue,
        m_maxNumRowsPerTrainingValue, m_numberOfTrainingValues, m_numberOfForecastValues;


    float m_minChangeAlpha, m_maxChangeAlpha;
    QString m_matrixOutputFile;
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
