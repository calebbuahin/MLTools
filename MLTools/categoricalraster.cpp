#include "stdafx.h"
#include "mrvm.h"
#include <QDebug>

CategoricalRaster::CategoricalRaster(MRVMItem::IOType iotype, const QString& name)
  :CategoricalMRVMItem(iotype, name), m_validCell(NULL) , m_columnCount(0) , m_driver(NULL)
{
  
}

CategoricalRaster::~CategoricalRaster()
{
  if(m_validCell)
      delete[] m_validCell;

  m_validCell = NULL;

  if(m_driver)
  delete m_driver;
  
  m_driver = NULL;
}

float* CategoricalRaster::trainingValues(int row)
{
  if(row < m_forecastValuesAsString.count())
    {
      return readDataFromRaster(m_forecastValuesAsString[row]);
    }

  return NULL;
}

void CategoricalRaster::setTrainingValuesAsString(const QList<QString> &trainingValues)
{
  CategoricalMRVMItem::setTrainingValuesAsString(trainingValues);
}

void CategoricalRaster::setTrainingValues(int row, float *&values)
{
  if(row < m_trainingValuesAsString.count())
    {
      QString filePath = m_trainingValuesAsString[row];
      writeDataToRaster(filePath, values);
    }
}

float* CategoricalRaster::forecastValues(int row)
{
  if(row < m_forecastValuesAsString.count())
    {
      return readDataFromRaster(m_forecastValuesAsString[row]);
    }

  return NULL;
}

void CategoricalRaster::setForecastValues(int row, float *&values)
{
  //check if file exists. otherwise create.
  if(row < m_forecastValuesAsString.count())
    {
      QString filePath = m_forecastValuesAsString[row];
      writeDataToRaster(filePath, values);
    }
}

void CategoricalRaster::setForecastValuesAsString(const QList<QString> &forecastValues)
{
  CategoricalMRVMItem::setForecastValuesAsString(forecastValues);
}

void CategoricalRaster::setForecastUncertaintyValueAsString(const QList<QString> &forecastUncertaintyValuesAsString)
{
    CategoricalMRVMItem::setForecastUncertaintyValueAsString(forecastUncertaintyValuesAsString);
}

float* CategoricalRaster::forecastUncertaintyValues(int row)
{

  if(row < m_forecastUncertaintyValuesAsString.count())
    {
      return readDataFromRaster(m_forecastUncertaintyValuesAsString[row]);
    }

  return NULL;
}

void CategoricalRaster::setForecastUncertaintyValues(int row, float *& values)
{
  //check if file exists. otherwise create.
  if(row < m_forecastUncertaintyValuesAsString.count())
    {
      QString filePath = m_forecastUncertaintyValuesAsString[row];
      writeDataToRaster(filePath, values);
    }

}


void CategoricalRaster::readXML(QXmlStreamReader &xmlReader)
{
  MRVMItem::readXML(xmlReader);
  readRasterProperties();

  if(m_iotype == MRVMItem::Output)
    {
      createOutputRasters();
    }
}

QString CategoricalRaster::type() const
{
  return "CategoricalRaster";
}

void CategoricalRaster::writeDataToRaster(const QString& filePath, float*& values)
{
  //check if file exists. otherwise create.

  GDALDataset* dataset = NULL;

  if(!QFile::exists(filePath) && m_iotype == MRVMItem::Output)
    {
      dataset = m_driver->Create(filePath.toStdString().c_str(), m_xSize , m_ySize, 1, GDT_Float32, NULL);
      dataset->SetGeoTransform(m_gcp);
      dataset->SetProjection(m_wktproj);
    }
  else
    {
      dataset = (GDALDataset*)GDALOpen(filePath.toStdString().c_str() , GA_Update);
    }

  if(dataset)
    {
      GDALRasterBand* dataBand = dataset->GetRasterBand(1);
      float * data = (float*) CPLMalloc(sizeof(float)*m_xSize*m_ySize);

      int count = 0;

      for(int i = 0 ; i < m_xSize ; i++)
        {
          for(int j = 0 ; j < m_ySize ; j++)
            {
              if(m_validCell[j * m_xSize + i])
                {
                  data[j*m_xSize + i] = values[count];
                  count++;
                }
              else
                {
                  data[j*m_xSize + i] = m_noData;
                }
            }
        }

      dataBand->RasterIO(GF_Write, 0,0,m_xSize , m_ySize, data, m_xSize , m_ySize , GDT_Float32, 0,0 );

      CPLFree(data);
      GDALClose(dataset);
      dataBand = NULL;
      dataset = NULL;
    }

  delete[] values;
}

float* CategoricalRaster::readDataFromRaster(const QString& filePath)
{
  float* values = NULL;

  if(QFile::exists(filePath) && m_columnCount)
    {
      GDALDataset* dataset = (GDALDataset*)GDALOpen(filePath.toStdString().c_str() , GA_ReadOnly);

      if(dataset)
        {
          GDALRasterBand* dataBand = dataset->GetRasterBand(1);
          float * data = (float*) CPLMalloc(sizeof(float)*m_xSize*m_ySize);
          values = new float[m_columnCount];

          dataBand->RasterIO(GF_Read, 0,0,m_xSize , m_ySize, data, m_xSize , m_ySize , GDT_Float32, 0,0 );


          int count = 0;

          for(int i = 0 ; i < m_xSize ; i++)
            {
              for(int j = 0 ; j < m_ySize ; j++)
                {
                  if(m_validCell[j * m_xSize + i])
                    {
                      values[count] = data[j*m_xSize + i];
                      count++;
                    }
                }
            }



          CPLFree(data);
          GDALClose(dataset);
          dataBand = NULL;
          dataset = NULL;

          return values;
        }
    }

  delete[] values;
}

void CategoricalRaster::readRasterProperties()
{
  if(m_trainingValuesAsString.count() > 0)
    {

      GDALDataset * dataset = (GDALDataset*)GDALOpen(m_trainingValuesAsString[0].toStdString().c_str(), GA_ReadOnly);
      m_driver =  dataset->GetDriver();
      dataset->GetGeoTransform(m_gcp);

      if(dataset)
        {
          assert(dataset->GetRasterCount() > 0);

          GDALRasterBand* rasterBand =  dataset->GetRasterBand(1);

          qDebug() << "Raster Type" << rasterBand->GetRasterDataType() ;

          m_xSize = dataset->GetRasterXSize();
          m_ySize = dataset->GetRasterYSize();
          m_noData = rasterBand->GetNoDataValue();
          m_wktproj = dataset->GetGCPProjection();

          float * data = (float*) CPLMalloc(sizeof(float)*m_xSize*m_ySize);

          rasterBand->RasterIO( GF_Read, 0,0,m_xSize , m_ySize, data, m_xSize , m_ySize , GDT_Float32, 0,0 );

          if(m_validCell)
            {
              delete[] m_validCell;
              m_validCell = NULL;
            }

          m_validCell = new int[m_xSize *  m_ySize];
          m_columnCount = 0;

          for(int i = 0 ; i < m_xSize ; i++)
            {
              for(int j = 0 ; j < m_ySize ; j++)
                {
                  if(data[j * m_xSize + i] == m_noData)
                    {
                      m_validCell[j*m_xSize + i] = 0;
                    }
                  else
                    {
                      m_validCell[j*m_xSize + i] = 1;
                      m_columnCount++;
                    }
                }
            }

          CPLFree(data);
          GDALClose(dataset);

          rasterBand = NULL;
          dataset = NULL;
        }
    }
}

void CategoricalRaster::createOutputRasters()
{
  if(m_driver)
    {
      for(int i = 0 ; i < m_forecastValuesAsString.count() ; i++)
        {
          GDALDataset* newData = m_driver->Create(m_forecastValuesAsString[0].toStdString().c_str() ,m_xSize , m_ySize , 1, GDT_CFloat32, NULL);
          newData->SetProjection(m_wktproj);
          newData->SetGeoTransform(m_gcp);
          GDALRasterBand* newBand = newData->GetRasterBand(1);
          newBand->SetNoDataValue(m_noData);
          GDALClose(newData);
        }

      for(int i = 0 ; i < m_forecastUncertaintyValuesAsString.count() ; i++)
        {
          GDALDataset* newData = m_driver->Create(m_forecastUncertaintyValuesAsString[0].toStdString().c_str() ,m_xSize , m_ySize , 1, GDT_CFloat32,NULL);
          newData->SetProjection(m_wktproj);
          newData->SetGeoTransform(m_gcp);
          GDALRasterBand* newBand = newData->GetRasterBand(1);
          newBand->SetNoDataValue(m_noData);
          GDALClose(newData);
        }
    }
}



