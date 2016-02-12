//#include <headers/stdafx.h>
//#include <headers/mrvm.h>
//#include <QDebug>

//CategoricalRaster::CategoricalRaster(MRVMItem::IOType iotype, const QString& name)
//  :CategoricalMRVMItem(iotype, name)
//{
//  m_validCell = NULL;
//  m_columnCount = 0;
//  m_numRowsPerTrainingValue = 0;
//  m_numRowsPerForecastValue = 0;
//  m_driver = NULL;
//}

//CategoricalRaster::~CategoricalRaster()
//{
//  if(m_validCell)
//    delete[] m_validCell;

//  m_validCell = NULL;

//  if(m_driver)
//    delete m_driver;

//  m_driver = NULL;
//}

//QString CategoricalRaster::getName()  const
//{
//  return m_name;
//}

//af::array CategoricalRaster::trainingValues(int row)
//{
//  if(row < m_trainingValuesAsString.count())
//    {
//      return readDataFromRaster(m_trainingValuesAsString[row]);
//    }

//  return af::array();
//}

//void CategoricalRaster::setTrainingValuesAsString(const QList<QString> &trainingValues)
//{
//  CategoricalMRVMItem::setTrainingValuesAsString(trainingValues);
//}

//af::array CategoricalRaster::forecastValues(int row)
//{
//  if(row < m_forecastValuesAsString.count())
//    {
//      return readDataFromRaster(m_forecastValuesAsString[row]);
//    }

//  return af::array();
//}

//void CategoricalRaster::setForecastValues(int row, const af::array& valuesf , const af::array& uncertf)
//{
//  //check if file exists. otherwise create.
//  if(row < m_forecastValuesAsString.count() && row < m_forecastUncertaintyValuesAsString.count() && MRVMItem::Output)
//    {
//      QString filePathForecast = m_forecastValuesAsString[row];
//      QString filePathForecastUncert = m_forecastUncertaintyValuesAsString[row];

//      GDALDataset* datasetForecast = NULL;
//      GDALDataset* datasetUncertainty = NULL;

//      if(!QFile::exists(filePathForecast))
//        {
//          datasetForecast = m_driver->Create(filePathForecast.toStdString().c_str(), m_xSize , m_ySize, 1, GDT_Float32, NULL);
//          datasetForecast->SetGeoTransform(m_gcp);
//          datasetForecast->SetProjection(m_wktproj);
//        }
//      else
//        {
//          datasetForecast = (GDALDataset*)GDALOpen(filePathForecast.toStdString().c_str() , GA_Update);
//        }

//      if(!QFile::exists(filePathForecastUncert))
//        {
//          datasetUncertainty = m_driver->Create(filePathForecastUncert.toStdString().c_str(), m_xSize , m_ySize, 1, GDT_Int32, NULL);
//          datasetUncertainty->SetGeoTransform(m_gcp);
//          datasetUncertainty->SetProjection(m_wktproj);
//        }
//      else
//        {
//          datasetUncertainty = (GDALDataset*)GDALOpen(filePathForecastUncert.toStdString().c_str() , GA_Update);
//        }


//      float* values = valuesf.host<float>();
//      float* uncert = uncertf.host<float>();


//      if(datasetForecast && datasetUncertainty)
//        {
//          int* classes = (int*) CPLMalloc(sizeof(int)*m_xSize*m_ySize);
//          float* classesuncert = (float*)  CPLMalloc(sizeof(float)*m_xSize*m_ySize);

//          int count = 0;

//          for(int i = 0 ; i < m_xSize ; i++)
//            {
//              for(int j = 0 ; j < m_ySize ; j++)
//                {
//                  if(m_validCell[j * m_xSize + i])
//                    {
//                      float maxValue = std::numeric_limits<float>::min();
//                      float uncertMaxValue = m_noData;

//                      int maxClass = 0;

//                      for(int c = 0 ; c < m_categorybyclass.count(); c++)
//                        {
//                          float tvalue =values[c *  m_numRowsPerValue + count];

//                          if(tvalue > maxValue)
//                            {
//                              maxValue = tvalue;
//                              uncertMaxValue = uncert[c *  m_numRowsPerValue + count];
//                              maxClass  = c;
//                            }
//                        }

//                      count++;
//                      int classIndex = m_classbyindex[maxClass];
//                      classes[j*m_xSize + i] = classIndex;
//                      classesuncert[j*m_xSize + i] = uncertMaxValue;
//                    }
//                  else
//                    {
//                      classes[j*m_xSize + i] = m_noData;
//                      classesuncert[j*m_xSize + i] = m_noData;
//                    }
//                }
//            }

//          GDALRasterBand* bandforecast = datasetForecast->GetRasterBand(1);
//          GDALRasterBand* bandforecastUncertainty = datasetUncertainty->GetRasterBand(1);

//          bandforecast->RasterIO(GF_Write, 0,0,m_xSize , m_ySize, classes, m_xSize , m_ySize , GDT_Int32, 0,0 );
//          bandforecastUncertainty->RasterIO(GF_Write, 0,0,m_xSize , m_ySize, classesuncert, m_xSize , m_ySize , GDT_Float32, 0,0 );

//          CPLFree(classes);
//          CPLFree(classesuncert);

//          GDALClose(datasetForecast);
//          GDALClose(datasetUncertainty);

//        }

//      delete[] values;
//      delete[] uncert;

//    }
//}

//void CategoricalRaster::setForecastValuesAsString(const QList<QString> &forecastValues)
//{
//  CategoricalMRVMItem::setForecastValuesAsString(forecastValues);
//}

//void CategoricalRaster::setForecastUncertaintyValueAsString(const QList<QString> &forecastUncertaintyValuesAsString)
//{
//  CategoricalMRVMItem::setForecastUncertaintyValueAsString(forecastUncertaintyValuesAsString);
//}

//void CategoricalRaster::readXML(QXmlStreamReader &xmlReader)
//{
//  CategoricalMRVMItem::readXML(xmlReader);
//  readRasterProperties();

//  if(m_iotype == MRVMItem::Output)
//    {
//      createOutputRasters();
//    }
//}

//int CategoricalRaster::numRowsPerTrainingValue() const
//{
//  return m_numRowsPerTrainingValue;
//}

//int CategoricalRaster::numRowsPerForecastValue() const
//{
//  return m_numRowsPerForecastValue;
//}

//int CategoricalRaster::columnCount()
//{
//  return m_categorybyclass.size();
//}

//QString CategoricalRaster::type() const
//{
//  return "CategoricalRaster";
//}

//af::array CategoricalRaster::readDataFromRaster(const QString& filePath)
//{
//  af::array values;

//  if(QFile::exists(filePath) && m_numRowsPerValue)
//    {
//      GDALDataset* dataset = (GDALDataset*)GDALOpen(filePath.toStdString().c_str() , GA_ReadOnly);

//      if(dataset)
//        {
//          GDALRasterBand* dataBand = dataset->GetRasterBand(1);
//          int * data = (int*) CPLMalloc(sizeof(int)*m_xSize*m_ySize);
//          values = af::array(m_numRowsPerValue, m_columnCount);

//          dataBand->RasterIO(GF_Read, 0,0,m_xSize , m_ySize, data, m_xSize , m_ySize , GDT_Int32, 0,0 );

//          int count = 0;

//          for(int i = 0 ; i < m_xSize ; i++)
//            {
//              for(int j = 0 ; j < m_ySize ; j++)
//                {
//                  if(m_validCell[j * m_xSize + i])
//                    {
//                      int pclass = data[j * m_xSize + i];

//                      for(int c = 0 ; c < m_categorybyclass.count(); c++)
//                        {
//                          if(m_classbyindex[c] == pclass)
//                            {
//                              values(count,c) = maxCValue;
//                            }
//                          else
//                            {
//                              values(count,c) = minCValue;
//                            }

//                        }

//                      count++;
//                    }
//                }
//            }

//          af_print(values);

//          CPLFree(data);
//          GDALClose(dataset);
//          dataBand = NULL;
//          dataset = NULL;

//          return values;
//        }
//    }

//  return values;
//}

//void CategoricalRaster::readRasterProperties()
//{
//  if(m_trainingValuesAsString.count() > 0)
//    {

//      GDALDataset * dataset = (GDALDataset*)GDALOpen(m_trainingValuesAsString[0].toStdString().c_str(), GA_ReadOnly);
//      m_driver =  dataset->GetDriver();
//      dataset->GetGeoTransform(m_gcp);

//      if(dataset)
//        {
//          ASSERT(dataset->GetRasterCount() > 0,"");

//          GDALRasterBand* rasterBand =  dataset->GetRasterBand(1);

//          qDebug() << "Raster Type" << rasterBand->GetRasterDataType() ;

//          m_xSize = dataset->GetRasterXSize();
//          m_ySize = dataset->GetRasterYSize();
//          m_noData = rasterBand->GetNoDataValue();
//          m_wktproj = dataset->GetGCPProjection();

//          int* data = (int*) CPLMalloc(sizeof(int)*m_xSize*m_ySize);

//          rasterBand->RasterIO(GF_Read, 0,0,m_xSize , m_ySize, data, m_xSize , m_ySize , GDT_Int32, 0,0 );

//          if(m_validCell)
//            {
//              delete[] m_validCell;
//              m_validCell = NULL;
//            }

//          m_validCell = new int[m_xSize *  m_ySize];
//          m_columnCount = m_classbycategory.size();
//          m_numRowsPerValue = 0;

//          for(int i = 0 ; i < m_xSize ; i++)
//            {
//              for(int j = 0 ; j < m_ySize ; j++)
//                {
//                  int pclass = data[j * m_xSize + i];

//                  if(pclass == (int)m_noData)
//                    {
//                      m_validCell[j*m_xSize + i] = 0;
//                    }
//                  else
//                    {
//                      QMap<int,int>::iterator it = m_indexbyclass.find(pclass);

//                      if(it != m_indexbyclass.end())
//                        {
//                          m_validCell[j*m_xSize + i] = 1;
//                          m_numRowsPerValue++;
//                        }
//                      else
//                        {

//                          m_validCell[j*m_xSize + i] = 0;
//                        }
//                    }
//                }
//            }

//          QVector<QPointF> bounds;

//          bounds.append(getCoordinates(QPoint(0,0)).toPoint());
//          bounds.append(getCoordinates(QPoint(m_xSize-1,0)).toPoint());
//          bounds.append(getCoordinates(QPoint(0,m_ySize -1)).toPoint());
//          bounds.append(getCoordinates(QPoint(m_xSize,m_ySize -1)).toPoint());

//          m_boundary = QPolygonF(bounds);

//          CPLFree(data);
//          GDALClose(dataset);

//          rasterBand = NULL;
//          dataset = NULL;
//        }
//    }
//}

//void CategoricalRaster::createOutputRasters()
//{
//  if(m_driver)
//    {
//      for(int i = 0 ; i < m_forecastValuesAsString.count() ; i++)
//        {
//          GDALDataset* newData = m_driver->Create(m_forecastValuesAsString[i].toStdString().c_str() ,m_xSize , m_ySize , 1, GDT_Int32, NULL);
//          newData->SetProjection(m_wktproj);
//          newData->SetGeoTransform(m_gcp);
//          GDALRasterBand* newBand = newData->GetRasterBand(1);
//          newBand->SetNoDataValue(m_noData);
//          GDALClose(newData);

//          newBand = NULL;
//          newData = NULL;

//        }

//      for(int i = 0 ; i < m_forecastUncertaintyValuesAsString.count() ; i++)
//        {
//          GDALDataset* newData = m_driver->Create(m_forecastUncertaintyValuesAsString[i].toStdString().c_str() ,m_xSize , m_ySize , 1, GDT_CFloat32,NULL);
//          newData->SetProjection(m_wktproj);
//          newData->SetGeoTransform(m_gcp);
//          GDALRasterBand* newBand = newData->GetRasterBand(1);
//          newBand->SetNoDataValue(m_noData);
//          GDALClose(newData);

//          newBand = NULL;
//          newData = NULL;
//        }
//    }
//}



