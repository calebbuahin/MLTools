#include "stdafx.h"
#include "mrvm.h"

CategoricalRaster::CategoricalRaster(const QString& name)
  :CategoricalMRVMItem(name), m_validCell(NULL) , m_columnCount(0) , m_driver(NULL)
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


double* CategoricalRaster::values(int index)
{
  if(index < m_values.count())
    {
      GDALDataset* dataset = (GDALDataset*)GDALOpen(m_values[index].toStdString().c_str(), GA_ReadOnly);

      if(dataset)
        {
          GDALRasterBand* rasterBand = dataset->GetRasterBand(1);
          int * data = (int*) CPLMalloc(sizeof(int)*m_xSize*m_ySize);

          rasterBand->RasterIO( GF_Read, 0,0,m_xSize , m_ySize, data, m_xSize , m_ySize , GDT_Int32, 0,0 );

          double* values = new double[m_columnCount];

          int count = 0;

          //GPU Candidate
          int catCount = m_categories.count();

          for(int i = 0 ; i < m_xSize ; i++)
            {
              for(int j = 0 ; j < m_ySize ; j++)
                {
                  if(m_validCell[j * m_xSize + i])
                    {
                      int dataV = data[j * m_xSize + i];

                      for(QMap<int,QString>::iterator k = m_categories.begin()
                          ; k != m_categories.end() ; k++)
                        {

                          if(k.key() == dataV)
                            {
                               values[count] = data[j*m_xSize + i] ;
                            }
                          else
                            {
                               values[count] = data[j*m_xSize + i] ;
                            }


                          count++;
                        }

                    }
                }
            }


          return values;

          CPLFree(data);
          GDALClose(dataset);

          rasterBand = NULL;
          dataset = NULL;
        }
    }

  return NULL;
}

void CategoricalRaster::addValue(const QString& value)
{
  MRVMItem::addValue(value);
  readRasterProperties();
}

void CategoricalRaster::setValues(int index, double *&values)
{
  //check if file exists. otherwise create.
  if(index < m_values.count())
    {
      QString filePath = m_values[index];

      GDALDataset* dataset = NULL;

      if(!QFile::exists(filePath))
        {
          dataset = m_driver->Create(filePath.toStdString().c_str(), m_xSize , m_ySize, 1, GDT_Int32, NULL);
          dataset->SetGeoTransform(m_gcp);
          dataset->SetProjection(m_wktproj.toStdString().c_str());
        }
      else
        {
          dataset = (GDALDataset*)GDALOpen(filePath.toStdString().c_str() , GA_Update);
        }

      if(dataset)
        {
          GDALRasterBand* dataBand = dataset->GetRasterBand(1);
          int * data = (int*) CPLMalloc(sizeof(int)*m_xSize*m_ySize);
          dataBand->RasterIO( GF_Read, 0,0,m_xSize , m_ySize, data, m_xSize , m_ySize , GDT_Int32, 0,0 );

          int count = 0;

          //GPU Candidate
          int catCount = m_categories.count();

          for(int i = 0 ; i < m_xSize ; i++)
            {
              for(int j = 0 ; j < m_ySize ; j++)
                {
                  if(m_validCell[j * m_xSize + i])
                    {
                      int dataV = m_noData;

                      double max = std::numeric_limits<double>::min();

                      for(QMap<int,QString>::iterator k = m_categories.begin()
                          ; k != m_categories.end() ; k++)
                        {
                          if(values[count] > max)
                            {
                              max = values[count];
                              dataV = k.key();
                            }

                          count++;
                        }

                      data[j * m_xSize + i] = dataV;
                    }
                  else
                    {
                      data[j * m_xSize + i] = m_noData;
                    }
                }
            }

          dataBand->RasterIO(GF_Write, 0,0,m_xSize , m_ySize, data, m_xSize , m_ySize , GDT_Int32, 0,0 );

          CPLFree(data);
          GDALClose(dataset);

          dataBand = NULL;
          dataset = NULL;
        }
    }
}

int CategoricalRaster::columnCount()
{
  if(!m_columnCount && m_values.length() > 0)
    {
      readRasterProperties();
    }

  return m_columnCount;
}

QString CategoricalRaster::type() const
{
    return "CategoricalRaster";
}

void CategoricalRaster::readRasterProperties()
{
  if(m_values.count() > 0)
    {
      QFileInfo file(m_values[0]);

      GDALDataset * dataset = (GDALDataset*)GDALOpen(file.absoluteFilePath().toStdString().c_str(), GA_ReadOnly);
      m_driver = dataset->GetDriver();
      dataset->GetGeoTransform(m_gcp);

      if(dataset)
        {
          assert(dataset->GetRasterCount() > 0);

          GDALRasterBand* rasterBand =  dataset->GetRasterBand(1);

          m_xSize = dataset->GetRasterXSize();
          m_ySize = dataset->GetRasterYSize();
          m_noData = rasterBand->GetNoDataValue();
          m_wktproj = QString(dataset->GetGCPProjection());

          int * data = (int*) CPLMalloc(sizeof(int)*m_xSize*m_ySize);

          rasterBand->RasterIO( GF_Read, 0,0,m_xSize , m_ySize, data, m_xSize , m_ySize , GDT_Int32, 0,0 );

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
                  int dataV = data[j * m_xSize + i];

                  if(dataV == m_noData || m_categories.find(dataV) == m_categories.end())
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

          m_columnCount = m_columnCount * m_categories.count();


          CPLFree(data);
          GDALClose(dataset);

          rasterBand = NULL; dataset = NULL;
        }
    }
}
