#include <include/stdafx.h>
#include <include/mrvm.h>
#include <iostream>
#include <QDebug>

CategoricalRaster::CategoricalRaster(MRVMItem::IOType iotype, const QString& name)
    :CategoricalMRVMItem(iotype, name) , RasterItem()
{
    m_validCell = NULL;

    m_properties["IncludeLocation"] = true;

    m_columnCount = 0;
    m_numRowsPerTrainingValue = 0;
    m_numRowsPerForecastValue = 0;
    m_driver = NULL;
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

QString CategoricalRaster::getName()  const
{
    return m_name;
}

af::array CategoricalRaster::trainingValues(int row)
{
    if(row < m_trainingValuesAsString.count())
    {
        if(m_useRasterBootstrap)
        {
            return readTrainingDataFromSampler(m_trainingValuesAsString[row]);
        }
        else
        {
            return readDataFromRaster(m_trainingValuesAsString[row]);
        }
    }

    return af::array();
}

void CategoricalRaster::setTrainingValuesAsString(const QList<QString> &trainingValues)
{
    CategoricalMRVMItem::setTrainingValuesAsString(trainingValues);
}

af::array CategoricalRaster::forecastValues(int row)
{
    if(row < m_forecastValuesAsString.count())
    {
        if(m_useRasterBootstrap)
        {
            return readForecastDataFromSampler(m_trainingValuesAsString[row]);
        }
        else
        {
            return readDataFromRaster(m_forecastValuesAsString[row]);
        }
    }

    return af::array();
}

void CategoricalRaster::setForecastValues(int row, const af::array& valuesf , const af::array& uncertf)
{
    //check if file exists. otherwise create.
    if(row < m_forecastValuesAsString.count() && row < m_forecastUncertaintyValuesAsString.count() && MRVMItem::Output)
    {
        QString filePathForecast = m_forecastValuesAsString[row];
        QString filePathForecastUncert = m_forecastUncertaintyValuesAsString[row];

        GDALDataset* datasetForecast = NULL;
        GDALDataset* datasetUncertainty = NULL;

        if(!QFile::exists(filePathForecast))
        {
            datasetForecast = m_driver->Create(filePathForecast.toStdString().c_str(), m_xSize , m_ySize, 1, GDT_Float32, NULL);
            datasetForecast->SetGeoTransform(m_gcp);
            datasetForecast->SetProjection(m_wktproj);
        }
        else
        {
            datasetForecast = (GDALDataset*)GDALOpen(filePathForecast.toStdString().c_str() , GA_Update);
        }

        if(!QFile::exists(filePathForecastUncert))
        {
            datasetUncertainty = m_driver->Create(filePathForecastUncert.toStdString().c_str(), m_xSize , m_ySize, 1, GDT_Int32, NULL);
            datasetUncertainty->SetGeoTransform(m_gcp);
            datasetUncertainty->SetProjection(m_wktproj);
        }
        else
        {
            datasetUncertainty = (GDALDataset*)GDALOpen(filePathForecastUncert.toStdString().c_str() , GA_Update);
        }


        if(datasetForecast && datasetUncertainty)
        {
            int* classes = (int*) CPLMalloc(sizeof(int)*m_xSize*m_ySize);
            float* classesuncert = (float*)  CPLMalloc(sizeof(float)*m_xSize*m_ySize);

            int count = 0;

            for(int i = 0 ; i < m_xSize ; i++)
            {
                for(int j = 0 ; j < m_ySize ; j++)
                {
                    if(m_validCell[j * m_xSize + i])
                    {
                        af::array max;
                        af::array indexes;
                        af::array tval = valuesf(count, af::span);
                        af::max(max,indexes,tval,1);
                        uint32_t* value = indexes.host<uint32_t>();
                        float* uncert = uncertf(count,value[0]).host<float>();
                        int classn = m_classbyindex[value[0]];
                        classes[j*m_xSize + i] = classn;
                        classesuncert[j*m_xSize + i] = uncert[0];

                        delete[] value;
                        delete[] uncert;

                        count++;
                    }
                    else
                    {
                        classes[j*m_xSize + i] = m_noData;
                        classesuncert[j*m_xSize + i] = m_noData;
                    }
                }
            }

            GDALRasterBand* bandforecast = datasetForecast->GetRasterBand(1);
            GDALRasterBand* bandforecastUncertainty = datasetUncertainty->GetRasterBand(1);

            bandforecast->RasterIO(GF_Write, 0,0,m_xSize , m_ySize, classes, m_xSize , m_ySize , GDT_Int32, 0,0 );
            bandforecastUncertainty->RasterIO(GF_Write, 0,0,m_xSize , m_ySize, classesuncert, m_xSize , m_ySize , GDT_Float32, 0,0 );

            CPLFree(classes);
            CPLFree(classesuncert);

            GDALClose(datasetForecast);
            GDALClose(datasetUncertainty);

        }
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

void CategoricalRaster::readXML(QXmlStreamReader &xmlReader)
{
    CategoricalMRVMItem::readXML(xmlReader);
    m_includeLocation =  m_properties["IncludeLocation"].toBool();

    readRasterProperties();

    if(m_iotype == MRVMItem::Output)
    {
        createOutputRasters();
    }
}

int CategoricalRaster::numRowsPerTrainingValue() const
{
    return m_numRowsPerTrainingValue;
}

int CategoricalRaster::numRowsPerForecastValue() const
{
    return m_numRowsPerForecastValue;
}

void CategoricalRaster::resetProperties()
{

    if(m_useRasterBootstrap)
    {
        m_numRowsPerTrainingValue = m_sampleLocations.length();
        m_numRowsPerForecastValue = m_numValidPixels;
    }
    else
    {
        m_numRowsPerTrainingValue = m_numValidPixels;
        m_numRowsPerForecastValue = m_numValidPixels;
    }

    if(m_includeLocation)
    {
        m_columnCount = m_categorybyclass.size() + 2;
    }
    else
    {
        m_columnCount = m_categorybyclass.size();
    }
}

int CategoricalRaster::columnCount()
{
    return m_columnCount;
}

QString CategoricalRaster::type() const
{
    return "CategoricalRaster";
}

af::array CategoricalRaster::readDataFromRaster(const QString& filePath)
{
    af::array values;

    if(QFile::exists(filePath) && m_numRowsPerTrainingValue)
    {
        GDALDataset* dataset = (GDALDataset*)GDALOpen(filePath.toStdString().c_str() , GA_ReadOnly);

        if(dataset)
        {
            GDALRasterBand* dataBand = dataset->GetRasterBand(1);
            int * data = (int*) CPLMalloc(sizeof(int)*m_xSize*m_ySize);
            values = af::constant(minCValue, m_numRowsPerTrainingValue, m_columnCount);

            dataBand->RasterIO(GF_Read, 0,0,m_xSize , m_ySize, data, m_xSize , m_ySize , GDT_Int32, 0,0 );

            if(m_includeLocation)
            {
                int count = 0;

                for(int i = 0 ; i < m_xSize ; i++)
                {
                    for(int j = 0 ; j < m_ySize ; j++)
                    {
                        if(m_validCell[j * m_xSize + i])
                        {
                            int pclass = data[j * m_xSize + i];
                            QPointF p = getCoordinates(i,j);
                            values(count,m_indexbyclass[pclass]) = maxCValue ;
                            values(count , m_categorybyclass.size()) =p.x();
                            values(count, m_categorybyclass.size() + 1) = p.y();

                            count++;
                        }
                    }
                }
            }
            else
            {
                int count = 0;

                for(int i = 0 ; i < m_xSize ; i++)
                {
                    for(int j = 0 ; j < m_ySize ; j++)
                    {
                        if(m_validCell[j * m_xSize + i])
                        {
                            int pclass = data[j * m_xSize + i];
                            values(count,m_indexbyclass[pclass]) = maxCValue ;
                            count++;
                        }
                    }
                }
            }

            //af_print(values);

            CPLFree(data);
            GDALClose(dataset);
            dataBand = NULL;
            dataset = NULL;

            return values;
        }
    }

    return values;
}

af::array CategoricalRaster::readTrainingDataFromSampler(const QString& filePath)
{
    af::array values = af::constant(minCValue, m_numRowsPerTrainingValue, m_columnCount);

    if(QFile::exists(filePath) && m_numRowsPerTrainingValue)
    {
        GDALDataset* dataset = (GDALDataset*)GDALOpen(filePath.toStdString().c_str() , GA_ReadOnly);

        if(dataset)
        {
            GDALRasterBand* dataBand = dataset->GetRasterBand(1);
            int * data = (int*) CPLMalloc(sizeof(int)*m_xSize*m_ySize);

            dataBand->RasterIO(GF_Read, 0,0,m_xSize , m_ySize, data, m_xSize , m_ySize , GDT_Int32, 0,0 );

            if(m_includeLocation)
            {
                for(int i = 0 ; i < m_numRowsPerTrainingValue ; i++)
                {
                    QPoint pp = m_sampleLocations[i];
                    int pclass = data[pp.y() * m_xSize + pp.x()];
                    QPointF p = getCoordinates(pp);
                    values(i,m_indexbyclass[pclass]) = maxCValue ;
                    values(i , m_categorybyclass.size()) = p.x();
                    values(i, m_categorybyclass.size() + 1) = p.y();

                }
            }
            else
            {
                for(int i = 0 ; i < m_numRowsPerTrainingValue ; i++)
                {
                    QPoint pp = m_sampleLocations[i];
                    int pclass = data[pp.y() * m_xSize + pp.x()];
                    values(i,m_indexbyclass[pclass]) = maxCValue ;

                }
            }

            //af_print(values);

            CPLFree(data);
            GDALClose(dataset);
            dataBand = NULL;
            dataset = NULL;

            return values;
        }
    }

    return values;
}

af::array CategoricalRaster::readForecastDataFromSampler(const QString& filePath)
{
    af::array values;

    if(QFile::exists(filePath) && m_numRowsPerForecastValue)
    {
        GDALDataset* dataset = (GDALDataset*)GDALOpen(filePath.toStdString().c_str() , GA_ReadOnly);

        if(dataset)
        {
            GDALRasterBand* dataBand = dataset->GetRasterBand(1);
            int * data = (int*) CPLMalloc(sizeof(int)*m_xSize*m_ySize);
            values = af::constant(minCValue, m_numRowsPerForecastValue, m_columnCount);
           
            
            dataBand->RasterIO(GF_Read, 0,0,m_xSize , m_ySize, data, m_xSize , m_ySize , GDT_Int32, 0,0 );

            if(m_includeLocation)
            {
                int count = 0;

                for(int i = 0 ; i < m_xSize ; i++)
                {
                    for(int j = 0 ; j < m_ySize ; j++)
                    {
                        if(m_validCell[j * m_xSize + i])
                        {
                            int pclass = data[j * m_xSize + i];
                            QPointF p = getCoordinates(i,j);
                            values(count , m_indexbyclass[pclass]) = maxCValue;
                            values(count , m_categorybyclass.size()) =p.x();
                            values(count, m_categorybyclass.size() + 1) = p.y();

                            count++;
                        }
                    }
                }
            }
            else
            {
                int count = 0;

                for(int i = 0 ; i < m_xSize ; i++)
                {
                    for(int j = 0 ; j < m_ySize ; j++)
                    {
                        if(m_validCell[j * m_xSize + i])
                        {
                            int pclass = data[j * m_xSize + i];
                            values(count , m_indexbyclass[pclass]) = maxCValue;
                            count++;
                        }
                    }
                }
            }

           // af_print(values);

            CPLFree(data);
            GDALClose(dataset);
            dataBand = NULL;
            dataset = NULL;

            return values;
        }
    }

    return values;
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
            ASSERT(dataset->GetRasterCount() > 0,"");

            GDALRasterBand* rasterBand =  dataset->GetRasterBand(1);

            qDebug() << "Raster Type" << rasterBand->GetRasterDataType() ;

            m_xSize = dataset->GetRasterXSize();
            m_ySize = dataset->GetRasterYSize();
            m_noData = rasterBand->GetNoDataValue();
            m_wktproj = dataset->GetGCPProjection();

            int* data = (int*) CPLMalloc(sizeof(int)*m_xSize*m_ySize);

            rasterBand->RasterIO(GF_Read, 0,0,m_xSize , m_ySize, data, m_xSize , m_ySize , GDT_Int32, 0,0 );

            if(m_validCell)
            {
                delete[] m_validCell;
                m_validCell = NULL;
            }

            m_validCell = new int[m_xSize *  m_ySize];
            m_columnCount = m_classbycategory.size();
            m_numValidPixels = 0;

            for(int i = 0 ; i < m_xSize ; i++)
            {
                for(int j = 0 ; j < m_ySize ; j++)
                {
                    int pclass = data[j * m_xSize + i];

                    if(pclass == (int)m_noData)
                    {
                        m_validCell[j*m_xSize + i] = 0;
                    }
                    else
                    {
                        QMap<int,int>::iterator it = m_indexbyclass.find(pclass);

                        if(it != m_indexbyclass.end())
                        {
                            m_validCell[j*m_xSize + i] = 1;
                            m_numValidPixels++;
                        }
                        else
                        {

                            m_validCell[j*m_xSize + i] = 0;
                        }
                    }
                }
            }

            QVector<QPointF> bounds;

            bounds.append(getCoordinates(QPoint(0,0)).toPoint());
            bounds.append(getCoordinates(QPoint(m_xSize-1,0)).toPoint());
            bounds.append(getCoordinates(QPoint(0,m_ySize -1)).toPoint());
            bounds.append(getCoordinates(QPoint(m_xSize,m_ySize -1)).toPoint());

            m_boundary = QPolygonF(bounds);

            CPLFree(data);
            GDALClose(dataset);

            rasterBand = NULL;
            dataset = NULL;

            resetProperties();
        }
    }
}

void CategoricalRaster::createOutputRasters()
{
    if(m_driver)
    {
        for(int i = 0 ; i < m_forecastValuesAsString.count() ; i++)
        {
            GDALDataset* newData = m_driver->Create(m_forecastValuesAsString[i].toStdString().c_str() ,m_xSize , m_ySize , 1, GDT_Int32, NULL);
            newData->SetProjection(m_wktproj);
            newData->SetGeoTransform(m_gcp);
            GDALRasterBand* newBand = newData->GetRasterBand(1);
            newBand->SetNoDataValue(m_noData);
            GDALClose(newData);

            newBand = NULL;
            newData = NULL;

        }

        for(int i = 0 ; i < m_forecastUncertaintyValuesAsString.count() ; i++)
        {
            GDALDataset* newData = m_driver->Create(m_forecastUncertaintyValuesAsString[i].toStdString().c_str() ,m_xSize , m_ySize , 1, GDT_CFloat32,NULL);
            newData->SetProjection(m_wktproj);
            newData->SetGeoTransform(m_gcp);
            GDALRasterBand* newBand = newData->GetRasterBand(1);
            newBand->SetNoDataValue(m_noData);
            GDALClose(newData);

            newBand = NULL;
            newData = NULL;
        }
    }
}



