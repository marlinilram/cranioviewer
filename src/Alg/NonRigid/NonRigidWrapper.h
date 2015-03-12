#ifndef NonRigidWrapper_H
#define NonRigidWrapper_H

#include <QObject>
#include "NonRigid.h"

class NonRigidWrapper : public QObject
{
    Q_OBJECT

public:
    NonRigidWrapper();
    ~NonRigidWrapper();

    NonRigid *getNonRigid() { return non_rigid; };
    void NonRigidIter();

signals:
    void updateRenderers();

private:
    NonRigid *non_rigid;
};

#endif