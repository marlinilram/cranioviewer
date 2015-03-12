#ifndef TrWidget_H
#define TrWidget_H

#include <iostream>
#include <QObject>
#include <QSlider>
#include <QDoubleSpinBox>

class TrWidget : public QObject
{
    Q_OBJECT

public:
    TrWidget();
    ~TrWidget();

    void setWidgets(QSlider *sliders[7], QDoubleSpinBox *spinBox[5]);
    void resetTrans();

private slots:
    void updateSliderTransLR(int value);
    void updateSpinBoxTransLR(double value);
    void updateSliderTransPA(int value);
    void updateSpinBoxTransPA(double value);
    void updateSliderTransIS(int value);
    void updateSpinBoxTransIS(double value);
    void updateTransMax(double value);
    void updateTransMin(double value);
    void updateSliderRotationLR(int value);
    void updateSliderRotationPA(int value);
    void updateSliderRotationIS(int value);
    void updateSliderScale(int value);
    void updateTr(int widget_id, double val);

signals:
    void updateTransform(double *cur_trans);

private:
    void updateSliderTrans(int value, int axial);
    void updateSpinBoxTrans(double value, int axial);

private:
    double tr_vals[7]; // Transform: X Y Z Rotate: X Y Z and Scale

    QSlider *slider_trans[3];
    QDoubleSpinBox *spinBox_trans[3];
    QSlider *slider_rot[3];
    QDoubleSpinBox *spinBox_trans_minmax[2];
    QSlider *slider_scal;
};

#endif