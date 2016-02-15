#include <include/stdafx.h>
#include <include/mrvm.h>
#include <random>

using namespace std;

RasterBootstrapSampler::RasterBootstrapSampler()
{
    
    
}

RasterBootstrapSampler::~RasterBootstrapSampler()
{
    
}

int RasterBootstrapSampler::numSamples() const
{
    return m_numSamples;
}

void RasterBootstrapSampler::setNumSamples(int numSamples)
{
    m_numSamples = numSamples;
}

QList<QString> RasterBootstrapSampler::rasterItemNames() const
{
    return m_rasterItemNames;
}

QMap<QString,RasterItem*> RasterBootstrapSampler::rasterItems() const
{
    return m_rasterItems;
}

void RasterBootstrapSampler::addRasterItem(RasterItem *rasterItem)
{
    m_rasterItems[rasterItem->getName()] = rasterItem;
    rasterItem->m_useRasterBootstrap = true;
}

bool RasterBootstrapSampler::removeRasterItem(RasterItem *rasteritem)
{
    rasteritem->m_useRasterBootstrap = false;
    return m_rasterItems.remove(rasteritem->getName());
}

QList<QPointF> RasterBootstrapSampler::samplingLocations() const
{
    return  m_samplingLocations;
}

void RasterBootstrapSampler::readXML(QXmlStreamReader &xmlReader)
{
    while (!(xmlReader.isEndElement() && !xmlReader.name().compare("RasterBootstrapSampler", Qt::CaseInsensitive)) && !xmlReader.hasError())
    {
        if(!xmlReader.name().compare("NumSamples" , Qt::CaseInsensitive))
        {
            m_numSamples = xmlReader.readElementText().toInt();
        }
        else if(!xmlReader.name().compare("RasterItems" , Qt::CaseInsensitive))
        {
            while (!(xmlReader.isEndElement() && !xmlReader.name().compare("RasterItems" , Qt::CaseInsensitive)) && !xmlReader.hasError())
            {
                if(!xmlReader.name().compare("RasterItem", Qt::CaseInsensitive))
                {
                    while (!(xmlReader.isEndElement() && !xmlReader.name().compare("RasterItem" , Qt::CaseInsensitive)) && !xmlReader.hasError())
                    {
                        QXmlStreamAttributes propAttributes =  xmlReader.attributes();

                        QString name;

                        for(QXmlStreamAttributes::iterator it = propAttributes.begin() ; it != propAttributes.end() ; it++)
                        {
                            if(!it->name().compare("Name", Qt::CaseInsensitive))
                            {
                                name = it->value().toString();
                            }
                         }

                        ASSERT(name.length(), "Invalid raster item name");
                        m_rasterItemNames.append(name);
                        m_rasterItems[name] = NULL;

                        xmlReader.readNext();
                    }
                }

                xmlReader.readNext() ;
            }
        }
        else if(!xmlReader.name().compare("SamplingLocations" , Qt::CaseInsensitive))
        {
            m_samplingLocations.clear();

            while (!(xmlReader.isEndElement() && !xmlReader.name().compare("SamplingLocations" , Qt::CaseInsensitive)) && !xmlReader.hasError())
            {
                if(!xmlReader.name().compare("Point" , Qt::CaseInsensitive))
                {
                    QPointF point;

                    while (!(xmlReader.isEndElement() && !xmlReader.name().compare("Point" , Qt::CaseInsensitive)) && !xmlReader.hasError())
                    {
                        if(!xmlReader.name().compare("X" , Qt::CaseInsensitive))
                        {
                            point.setX(xmlReader.readElementText().toFloat());
                        }
                        else if(!xmlReader.name().compare("Y" , Qt::CaseInsensitive))
                        {
                            point.setY(xmlReader.readElementText().toFloat());
                        }

                        xmlReader.readNext();
                    }

                    m_samplingLocations.append(point);
                }

                xmlReader.readNext();
            }
        }

        xmlReader.readNext();
    }
    
}

void RasterBootstrapSampler::writeXML(QXmlStreamWriter &xmlWriter)
{
    xmlWriter.writeStartElement("RasterBootstrapSampler");

    xmlWriter.writeTextElement("NumSamples", QString::number(m_numSamples));

    if(m_rasterItemNames.length())
    {
        xmlWriter.writeStartElement("RasterItems");

        for(int i = 0 ; i < m_rasterItemNames.length() ; i++)
        {
            QString name = m_rasterItemNames[i];

            xmlWriter.writeStartElement("RasterItem");

            xmlWriter.writeAttribute("Name", name);

            xmlWriter.writeEndElement();
        }

        xmlWriter.writeEndElement();
    }

    if(m_samplingLocations.count())
    {
        xmlWriter.writeStartElement("SamplingLocations");

        for(int i = 0 ; i < m_samplingLocations.length() ; i++)
        {
            QPointF point = m_samplingLocations[i];
            xmlWriter.writeStartElement("Point");
            xmlWriter.writeTextElement("X" , QString::number(point.x()));
            xmlWriter.writeTextElement("Y" , QString::number(point.y()));
            xmlWriter.writeEndElement();

        }

        xmlWriter.writeEndElement();
    }

    xmlWriter.writeEndElement();
}

void RasterBootstrapSampler::createValidWindowCenters()
{
    m_samplingLocations.clear();

    if(m_rasterItems.size())
    {

        //Overkill but ah well!
        RasterItem* item = m_rasterItems.first();

        QPolygonF p = item->boundary();

        for(QMap<QString,RasterItem*>::Iterator it = m_rasterItems.begin();
            it !=  m_rasterItems.end() ; it++)
        {
            if(it.value() != item)
            {
                RasterItem* item = it.value();
                p = p.intersected(item->boundary());
            }
        }

        QRectF rect = p.boundingRect();
        uniform_real_distribution<float> xdist(rect.left() , rect.right());
        uniform_real_distribution<float> ydist(rect.bottom() , rect.top());

        while (m_samplingLocations.length() < m_numSamples)
        {
            QPointF p;
            p.setX(xdist(gen));
            p.setY(ydist(gen));

            bool noerror = true ;

            for(QMap<QString,RasterItem*>::Iterator it = m_rasterItems.begin();
                it !=  m_rasterItems.end() ; it++)
            {
                RasterItem* item = it.value();
                QPoint index;

                if(!(item->contains(p,index) && item->isValid(index)))
                {
                    noerror = false;
                    break;
                }
            }

            if(noerror)
            {
                m_samplingLocations.append(p);
            }
        }
    }
}

void RasterBootstrapSampler::setRasterItemSamplingAttributes()
{
    for(QMap<QString,RasterItem*>::Iterator it = m_rasterItems.begin();
        it !=  m_rasterItems.end() ; it++)
    {
        RasterItem* item = it.value();
        item->setBootstrapSamplingPoints(m_samplingLocations);
    }
}
