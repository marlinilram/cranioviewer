#include "TrWidget.h"

TrWidget::TrWidget()
{
    for (size_t i = 0; i < 6; ++i)
    {
        tr_vals[i] = 0.0;
    }
    tr_vals[6] = 1.0;
}

TrWidget::~TrWidget()
{

}

void TrWidget::updateSliderTransLR(int value)
{
    updateSliderTrans(value, 0);
}

void TrWidget::updateSpinBoxTransLR(double value)
{
    updateSpinBoxTrans(value, 0);
}

void TrWidget::updateSliderTransPA(int value)
{
    updateSliderTrans(value, 1);
}

void TrWidget::updateSpinBoxTransPA(double value)
{
    updateSpinBoxTrans(value, 1);
}

void TrWidget::updateSliderTransIS(int value)
{
    updateSliderTrans(value, 2);
}

void TrWidget::updateSpinBoxTransIS(double value)
{
    updateSpinBoxTrans(value, 2);
}

void TrWidget::updateSliderTrans(int value, int axial)
{
    double new_val = spinBox_trans[axial]->minimum()+(spinBox_trans[axial]->maximum()-spinBox_trans[axial]->minimum())*value/(slider_trans[axial]->maximum()-slider_trans[axial]->minimum());
    spinBox_trans[axial]->setValue(new_val);
    updateTr(axial, new_val);
}

void TrWidget::updateSpinBoxTrans(double value, int axial)
{
    double new_val = (slider_trans[axial]->maximum()-slider_trans[axial]->minimum())*(value-spinBox_trans[axial]->minimum())/(spinBox_trans[axial]->maximum()-spinBox_trans[axial]->minimum())+slider_trans[axial]->minimum();
    slider_trans[axial]->setValue((slider_trans[axial]->maximum()-slider_trans[axial]->minimum())*(value-spinBox_trans[axial]->minimum())/(spinBox_trans[axial]->maximum()-spinBox_trans[axial]->minimum())+slider_trans[axial]->minimum());
    //updateTr(axial, new_val);
}

void TrWidget::updateTransMax(double value)
{
    for (size_t i = 0; i < 3; ++i)
    {
        slider_trans[i]->setRange(0, 100*(value-spinBox_trans_minmax[0]->value()));
        spinBox_trans[i]->setRange(spinBox_trans_minmax[0]->value(), spinBox_trans_minmax[1]->value());
    }
}

void TrWidget::updateTransMin(double value)
{
    for (size_t i = 0; i < 3; ++i)
    {
        slider_trans[i]->setRange(0, 100*(spinBox_trans_minmax[1]->value()-value));
        spinBox_trans[i]->setRange(spinBox_trans_minmax[0]->value(), spinBox_trans_minmax[1]->value());
    }
}

void TrWidget::updateSliderRotationLR(int value)
{
    updateTr(3, slider_rot[0]->value());
}

void TrWidget::updateSliderRotationPA(int value)
{
    updateTr(4, slider_rot[1]->value());
}

void TrWidget::updateSliderRotationIS(int value)
{
    updateTr(5, slider_rot[2]->value());
}

void TrWidget::updateSliderScale(int value)
{
    updateTr(6, (double)slider_scal->value()/100);
}

void TrWidget::updateTr(int widget_id, double val)
{
    tr_vals[widget_id] = val;

    emit(updateTransform(tr_vals));
}

void TrWidget::resetTrans()
{
    for (size_t i = 0; i < 6; ++i)
    {
        tr_vals[i] = 0.0;
    }
    tr_vals[6] = 1.0;
    emit(updateTransform(tr_vals));
}

void TrWidget::setWidgets(QSlider *sliders[7], QDoubleSpinBox *spinBox[5])
{
    for (size_t i = 0; i < 3; ++i)
    {
        slider_trans[i] = sliders[i];
        spinBox_trans[i] = spinBox[i];
    }
    for (size_t i = 3; i < 5; ++i)
    {
        spinBox_trans_minmax[i-3] = spinBox[i];
    }
    for (size_t i = 3; i < 6; ++i)
    {
        slider_rot[i-3] = sliders[i];
    }
    slider_scal = sliders[6];
}
