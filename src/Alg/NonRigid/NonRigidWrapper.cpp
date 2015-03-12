#include "NonRigidWrapper.h"

NonRigidWrapper::NonRigidWrapper()
{
    non_rigid = new NonRigid;
}

NonRigidWrapper::~NonRigidWrapper()
{
    delete non_rigid;
}
