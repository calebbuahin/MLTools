#include "mrvm.h"

RasterItem::RasterItem()
{

}

RasterItem::~RasterItem()
{

}

bool RasterItem::contains(const QPointF &point)
{
  QPoint f  = getCoordinateIndexes(point);

  if(f.x() < 0 || f.x() >= m_xSize || f.y() < 0 || f.y() >= m_ySize)
    return false;

  return true;
}

bool RasterItem::isValid(const QPoint& index)
{
  return m_validCell[index.y() * m_xSize + index.x()];
}

QPointF RasterItem::getCoordinates(const QPoint& indexes) const
{
  QPointF p;

  p.setX( m_gcp[0] + indexes.x() * m_gcp[1] + indexes.y() * m_gcp[2]);
  p.setY(m_gcp[3] + indexes.x() * m_gcp[4] + indexes.y() * m_gcp[5]);

  return p;
}

QPoint RasterItem::getCoordinateIndexes(const QPointF &coordinates) const
{
  QPoint p;

  p.setY((coordinates.y() - m_gcp[3] - ((coordinates.x() - m_gcp[0]) * m_gcp[4] / m_gcp[1]))/(m_gcp[5] - ( m_gcp[2] * m_gcp[4] / m_gcp[1])));
  p.setX((coordinates.x() - m_gcp[0]) / m_gcp[1] - (p.y() * m_gcp[2]) / m_gcp[1]);

  return p;
}

void RasterItem::setBootStrapPoints(const QList<QPoint> &centers, const QList<QList<QPoint> > &indexes)
{
  m_bootStrapCenters = centers;
  m_bootStrapSamplingPoints = indexes;
}

QPolygonF RasterItem::boundary() const
{
  QVector<QPointF> bounds;


  bounds.append(getCoordinates(QPoint(0,0)).toPoint());
  bounds.append(getCoordinates(QPoint(m_xSize-1,0)).toPoint());
  bounds.append(getCoordinates(QPoint(0,m_ySize -1)).toPoint());
  bounds.append(getCoordinates(QPoint(m_xSize,m_ySize -1)).toPoint());

  return QPolygonF(bounds);
}
