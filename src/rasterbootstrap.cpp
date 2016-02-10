#include "mrvm.h"
#include <random>

RasterBootstrap::RasterBootstrap()
{


}


RasterBootstrap::~RasterBootstrap()
{

}

int RasterBootstrap::numSamplingWindows() const
{
  return m_numSampleWindows;
}

void RasterBootstrap::setNumSamplingWindows(int windowCount)
{
  m_numSampleWindows = windowCount;
}

int RasterBootstrap::samplingWindowSize() const
{
  return m_samplingWindowSize;
}

void RasterBootstrap::setSamplingWindowSize(int size)
{
  m_samplingWindowSize = size;
}


bool RasterBootstrap::includeDistance() const
{
  return m_includeDistance;
}

void RasterBootstrap::setIncludeDistance(bool includedistance)
{
  m_includeDistance = includedistance;
}

QSet<RasterItem*> RasterBootstrap::rasterItems() const
{
  return m_rasterItems;
}

void RasterBootstrap::addRasterItem(RasterItem *rasterItem)
{
  m_rasterItems.insert(rasterItem);
}

bool RasterBootstrap::removeRasterItem(RasterItem *rasteriterm)
{
  return m_rasterItems.remove(rasteriterm);
}

QMap<QString,QList<QList<QPoint>>> RasterBootstrap::sampleLocationIndexes() const
{
  return  m_mrvmItemLocations;
}

QMap<QString,QList<QPoint>> RasterBootstrap::windowLocations() const
{
  return m_windowLocations;
}

void RasterBootstrap::sampleRasters()
{
  m_mrvmItemLocations.clear();
  m_windowLocations.clear();

  if(m_rasterItems.count())
    {

      QPolygonF poly = (*m_rasterItems.begin())->boundary();

      for(QSet<RasterItem*>::iterator it ; it != m_rasterItems.end() ; it++)
        {
          RasterItem* item = *it;
          poly = poly.intersected(item->boundary());
        }

      QRectF rect = poly.boundingRect();

      std::default_random_engine generator;
      std::uniform_real_distribution<double> distributionX(rect.left(),rect.right());
      std::uniform_real_distribution<double> distributionY(rect.top(),rect.bottom());

      while (m_windowLocations.size() < m_numSampleWindows)
        {
          double x = distributionX(generator);
          double y = distributionY(generator);

          if(poly.contains(QPoint(x,y)))
          {



          }
        }
    }
}

void RasterBootstrap::setRasterItemLocations()
{
  for(QSet<RasterItem*>::iterator it = m_rasterItems.begin() ;
      it != m_rasterItems.end() ; it++)
    {
      RasterItem* rasterItem = *it;

      if(m_mrvmItemLocations.contains(rasterItem->getName()) &&
         m_windowLocations.contains(rasterItem->getName()))
        {
          QList<QList<QPoint>> mrvmitems = m_mrvmItemLocations[rasterItem->getName()];
          QList<QPoint> windowLocations = m_windowLocations[rasterItem->getName()];
          rasterItem->setBootStrapPoints(windowLocations , mrvmitems );
        }
    }
}
