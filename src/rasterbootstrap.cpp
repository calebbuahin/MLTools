#include <headers/stdafx.h>
#include <headers/mrvm.h>
#include <random>

RasterBootstrapSampler::RasterBootstrapSampler()
{


}


RasterBootstrapSampler::~RasterBootstrapSampler()
{

}

int RasterBootstrapSampler::numSamplingWindows() const
{
  return m_numSampleWindows;
}

void RasterBootstrapSampler::setNumSamplingWindows(int windowCount)
{
  m_numSampleWindows = windowCount;
}

int RasterBootstrapSampler::samplingWindowSize() const
{
  return m_samplingWindowSize;
}

void RasterBootstrapSampler::setSamplingWindowSize(int size)
{
  m_samplingWindowSize = size;
}


bool RasterBootstrapSampler::includeDistance() const
{
  return m_includeDistance;
}

void RasterBootstrapSampler::setIncludeDistance(bool includedistance)
{
  m_includeDistance = includedistance;
}

QList<QString> RasterBootstrapSampler::rasterItemNames() const
{
    return m_rasterItemNames;
}

QMap<QString,int> RasterBootstrapSampler::sampleSizeForRasterItem() const
{
    return m_sampleSizeForRasterItems;
}

QMap<QString,RasterItem*> RasterBootstrapSampler::rasterItems() const
{
  return m_rasterItems;
}

void RasterBootstrapSampler::addRasterItem(RasterItem *rasterItem)
{
  m_rasterItems[rasterItem->getName()] = rasterItem;

}

bool RasterBootstrapSampler::removeRasterItem(RasterItem *rasteritem)
{
  rasteritem->m_useRasterBootstrap = false;

  return m_rasterItems.remove(rasteriterm->getName());
}

QMap<QString,QList<QList<QPointF>>> RasterBootstrapSampler::sampleLocationIndexes() const
{
  return  m_mrvmItemLocations;
}

QMap<QString,QList<QPointF>> RasterBootstrapSampler::windowLocations() const
{
  return m_windowLocations;
}

void RasterBootstrapSampler::sampleRasters()
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


