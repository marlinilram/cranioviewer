#ifndef DistMapConfigDialog_H
#define DistMapConfigDialog_H

#include <QDialog>
#include "ui_DistMapConfigDialog.h"

class DistMapConfigDialog : public QDialog, public Ui::DistMapConfigDialog
{
    Q_OBJECT

public:
    DistMapConfigDialog();
    ~DistMapConfigDialog();
};

#endif