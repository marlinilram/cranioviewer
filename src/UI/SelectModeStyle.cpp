#include "SelectModeStyle.h"

#include <vtkCommand.h>
#include <vtkNew.h>
#include <vtkObjectFactory.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkUnsignedCharArray.h>

vtkStandardNewMacro(SelectModeStyle);

SelectModeStyle::SelectModeStyle()
{
  StartPtId = -1;
}

SelectModeStyle::~SelectModeStyle()
{

}

void SelectModeStyle::OnMouseMove()
{
  if (!this->Interactor || !this->Moving)
  {
    return;
  }

  // pressed and move
  int NewEndPosition[2];
  NewEndPosition[0] = this->Interactor->GetEventPosition()[0];
  NewEndPosition[1] = this->Interactor->GetEventPosition()[1];
  int *size = this->Interactor->GetRenderWindow()->GetSize();
  if (NewEndPosition[0] > (size[0]-1))
  {
    NewEndPosition[0] = size[0]-1;
  }
  if (NewEndPosition[0] < 0)
  {
    NewEndPosition[0] = 0;
  }
  if (NewEndPosition[1] > (size[1]-1))
  {
    NewEndPosition[1] = size[1]-1;
  }
  if (NewEndPosition[1] < 0)
  {
    NewEndPosition[1] = 0;
  }

  if (sqrt((NewEndPosition[0]-this->EndPosition[0])*(NewEndPosition[0]-this->EndPosition[0])
    + (NewEndPosition[1]-this->EndPosition[1])*(NewEndPosition[1]-this->EndPosition[1])) > 10)
  {
    this->EndPosition[0] = NewEndPosition[0];
    this->EndPosition[1] = NewEndPosition[1];
    if(this->DrawPolygonPixels)
    {
      this->DrawPolygon();

    }
  }
}

void SelectModeStyle::OnLeftButtonDown()
{
  if (!this->Interactor)
  {
    return;
  }

  this->Moving = 1;

  vtkRenderWindow *renWin = this->Interactor->GetRenderWindow();

  this->StartPosition[0] = this->Interactor->GetEventPosition()[0];
  this->StartPosition[1] = this->Interactor->GetEventPosition()[1];
  this->EndPosition[0] = this->StartPosition[0];
  this->EndPosition[1] = this->StartPosition[1];

  this->PixelArray->Initialize();
  this->PixelArray->SetNumberOfComponents(3);
  int *size = renWin->GetSize();
  this->PixelArray->SetNumberOfTuples(size[0]*size[1]);

  renWin->GetPixelData(0, 0, size[0]-1, size[1]-1, 1, this->PixelArray);

  this->InvokeEvent(vtkCommand::StartInteractionEvent);
}

void SelectModeStyle::OnLeftButtonUp()
{
  if (!this->Interactor || !this->Moving)
  {
    return;
  }

  if(this->DrawPolygonPixels)
  {
    int *size = this->Interactor->GetRenderWindow()->GetSize();
    unsigned char *pixels = this->PixelArray->GetPointer(0);
    this->Interactor->GetRenderWindow()->SetPixelData(
      0, 0, size[0]-1, size[1]-1, pixels, 1);
  }

  this->Moving = 0;
  this->StartPtId = -1;
  this->InvokeEvent(vtkCommand::SelectionChangedEvent);
  this->InvokeEvent(vtkCommand::EndInteractionEvent);
}

void SelectModeStyle::DrawPolygon()
{
  vtkNew<vtkUnsignedCharArray> tmpPixelArray;
  tmpPixelArray->DeepCopy(this->PixelArray);
  unsigned char *pixels = tmpPixelArray->GetPointer(0);
  int *size = this->Interactor->GetRenderWindow()->GetSize();

  // draw a line from the end to the start
  this->DrawPixels(StartPosition, EndPosition, pixels, size);

  this->Interactor->GetRenderWindow()->SetPixelData(0, 0, size[0]-1, size[1]-1, pixels, 1);
}

void SelectModeStyle::DrawPixels(int StartPos[2], int EndPos[2], unsigned char *pixels, int *size)
{
  int x1=StartPos[0], x2=EndPos[0];
  int y1=StartPos[1], y2=EndPos[1];

  double x = x2 - x1;
  double y = y2 - y1;
  double length = sqrt( x*x + y*y );
  if(length == 0)
  {
    return;
  }
  double addx = x / length;
  double addy = y / length;

  x = x1;
  y = y1;
  int row, col;
  for(double i = 0; i < length; i += 1)
  {
    col = (int)x;
    row = (int)y;
    pixels[3*(row*size[0]+col)] = 255 ^ pixels[3*(row*size[0]+col)];
    pixels[3*(row*size[0]+col)+1] = 255 ^ pixels[3*(row*size[0]+col)+1];
    pixels[3*(row*size[0]+col)+2] = 255 ^ pixels[3*(row*size[0]+col)+2];
    x += addx;
    y += addy;
  }
}

void SelectModeStyle::GetStartPos(int StartPos[2])
{
  StartPos[0] = this->StartPosition[0];
  StartPos[1] = this->StartPosition[1];
}

void SelectModeStyle::GetEndPos(int EndPos[2])
{
  EndPos[0] = this->EndPosition[0];
  EndPos[1] = this->EndPosition[1];
}

double SelectModeStyle::GetRadius()
{
  return sqrt((this->StartPosition[0]-this->EndPosition[0])*(this->StartPosition[0]-this->EndPosition[0])
    + (this->StartPosition[1]-this->EndPosition[1])*(this->StartPosition[1]-this->EndPosition[1]));
}

int SelectModeStyle::GetMovingTag()
{
  return this->Moving;
}

void SelectModeStyle::SetStartPtId(vtkIdType id)
{
  this->StartPtId = id;
}

vtkIdType SelectModeStyle::GetStartPtId()
{
  return this->StartPtId;
}