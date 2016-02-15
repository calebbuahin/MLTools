#include <include/stdafx.h>
#include <include/mrvm.h>
#include <random>

using namespace std;

RasterItem::RasterItem()
{

}

RasterItem::~RasterItem()
{

}

bool RasterItem::contains(const QPointF &point, QPoint& pointIndex)
{
    pointIndex  = getCoordinateIndexes(point);

    if(pointIndex.x() < 0 || pointIndex.x() >= m_xSize || pointIndex.y() < 0 || pointIndex.y() >= m_ySize)
        return false;

    return true;
}

bool RasterItem::isValid(const QPoint& index)
{
    return m_validCell[index.y() * m_xSize + index.x()];
}

bool RasterItem::isValid(int x , int y)
{
    return m_validCell[y* m_xSize + x];
}

QPointF RasterItem::getCoordinates(const QPoint& indexes) const
{
    QPointF p;

    p.setX(m_gcp[0] + indexes.x() * m_gcp[1] + indexes.y() * m_gcp[2]);
    p.setY(m_gcp[3] + indexes.x() * m_gcp[4] + indexes.y() * m_gcp[5]);

    return p;
}

QPointF RasterItem::getCoordinates(int x, int y) const
{
    QPointF p;

    p.setX(m_gcp[0] + x * m_gcp[1] + y * m_gcp[2]);
    p.setY(m_gcp[3] + x * m_gcp[4] + y * m_gcp[5]);

    return p;
}

QPoint RasterItem::getCoordinateIndexes(const QPointF &coordinates) const
{
    QPoint p;

    p.setY((coordinates.y() - m_gcp[3] - ((coordinates.x() - m_gcp[0]) * m_gcp[4] / m_gcp[1]))/(m_gcp[5] - ( m_gcp[2] * m_gcp[4] / m_gcp[1])));
    p.setX((coordinates.x() - m_gcp[0]) / m_gcp[1] - (p.y() * m_gcp[2]) / m_gcp[1]);

    return p;
}

void RasterItem::setBootstrapSamplingPoints(const QList<QPointF>& sampleLocations)
{
    m_sampleLocations.clear();

    for(int i = 0 ; i < sampleLocations.length() ; i++)
    {
        QPoint p = getCoordinateIndexes(sampleLocations[i]);
        m_sampleLocations.append(p);
    }

    resetProperties();
}

QPolygonF RasterItem::boundary() const
{
    return m_boundary;
}

