//#include "stdafx.h"
//#include "mrvm.h"

//CategoricalRaster::CategoricalRaster(MRVMItem::IOType iotype, const QString& name)
//  :CategoricalMRVMItem(iotype, name), m_validCell(NULL) , m_columnCount(0) , m_driver(NULL)
//{
  
//}

//CategoricalRaster::~CategoricalRaster()
//{
//  if(m_validCell)
//      delete[] m_validCell;

//  m_validCell = NULL;

//  if(m_driver)
//  delete m_driver;
  
//  m_driver = NULL;
//}

//float* CategoricalRaster::trainingValues(int row)
//{
//  if(row < m_forecastValuesAsString.count())
//    {
//      return readDataFromRaster(m_forecastValuesAsString[row]);
//    }

//  return NULL;
//}

//void CategoricalRaster::setTrainingValuesAsString(const QList<QString> &trainingValues)
//{
//  CategoricalMRVMItem::setTrainingValuesAsString(trainingValues);
//}

//void CategoricalRaster::setTrainingValues(int row, float *&values)
//{
//  if(row < m_trainingValuesAsString.count())
//    {
//      QString filePath = m_trainingValuesAsString[row];
//      writeDataToRaster(filePath, values);
//    }
//}

//float* CategoricalRaster::forecastValues(int row)
//{
//  if(row < m_forecastValuesAsString.count())
//    {
//      return readDataFromRaster(m_forecastValuesAsString[row]);
//    }

//  return NULL;
//}

//void CategoricalRaster::setForecastValues(int row, float *&values)
//{
//  //check if file exists. otherwise create.
//  if(row < m_forecastValuesAsString.count())
//    {
//      QString filePath = m_forecastValuesAsString[row];
//      writeDataToRaster(filePath, values);
//    }
//}

//void CategoricalRaster::setForecastValuesAsString(const QList<QString> &forecastValues)
//{
//  CategoricalMRVMItem::setForecastValuesAsString(forecastValues);
//}

//void CategoricalRaster::setForecastUncertaintyValueAsString(const QList<QString> &forecastUncertaintyValuesAsString)
//{
//    CategoricalMRVMItem::setForecastUncertaintyValueAsString(forecastUncertaintyValuesAsString);
//}

//float* CategoricalRaster::forecastUncertaintyValues(int row)
//{

//  if(row < m_forecastUncertaintyValuesAsString.count())
//    {
//      return readDataFromRaster(m_forecastUncertaintyValuesAsString[row]);
//    }

//  return NULL;
//}

