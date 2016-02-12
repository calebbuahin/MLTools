#include <include/stdafx.h>
#include <include/mrvm.h>

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

void RasterItem::setBootstrapPoints(const QList<QPoint> & windowCenters, int samplingWindowSize, int numSamples, bool includeDistanceWithBootstrap)
{

}

bool RasterItem::includeDistanceWithBootstrap() const
{
    return m_includeDistanceWithBootstrap;
}

QPolygonF RasterItem::boundary() const
{
    return m_boundary;
}


QList<QPoint> RasterItem::sampleRasterForPointsWithinWindow(const QPoint &center)
{
    QList<QPoint> points
}
