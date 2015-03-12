#include "CranioViewer.h"
#include <QTextCodec>

int main( int argc, char **argv )
{
	QApplication *app = new QApplication(argc, argv);
//	QTextCodec::setCodecForTr(QTextCodec::codecForName("GB2312"));

	CranioViewer *window = new CranioViewer();
	window->show();
	return app->exec();
}