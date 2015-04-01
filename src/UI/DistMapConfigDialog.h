#ifndef DistMapConfigDialog_H
#define DistMapConfigDialog_H

#include <QDialog>
#include "ui_DistMapConfig.h"

class DistMapConfigDialog : public QDialog, public Ui::DistMapConfigDialog
{
    Q_OBJECT

public:
    DistMapConfigDialog();
    ~DistMapConfigDialog();
};

#endif