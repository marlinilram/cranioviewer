CMAKE_MINIMUM_REQUIRED( VERSION 2.8 )
PROJECT( CranioViewer )

SET( CMAKE_PREFIX_PATH "C:\\Qt\\Qt5.5.1\\5.5\\msvc2013_64" )


# QT FILES						
FIND_PACKAGE( Qt5Widgets )
FIND_PACKAGE( Qt5Xml )
#FIND_PACKAGE( Qt5OpenGL )
FIND_PACKAGE( Qt5Core )

# VTK FILES
FIND_PACKAGE( VTK REQUIRED )
INCLUDE( ${VTK_USE_FILE} )
SET( vtk_lib_renew ${VTK_LIBRARIES} )
#SET( vtk_use_renew ${VTK_USE_FILE} )

# ITK FILES
#FIND_PACKAGE( ITK REQUIRED )
#INCLUDE( ${ITK_USE_FILE} )

# BOOST FILES
FIND_PATH( BOOST_DIR "boost" )

# OTHER DEPENDENCIES
#FIND_PACKAGE( OpenGL REQUIRED )
#FIND_PACKAGE( GLEW REQUIRED )
#FIND_PACKAGE( OpenCV REQUIRED )

# SET PROJECT VARS
# SET SRC VARS
SET( Project_INCLUDE_DIR 	src/
							src/UI
							src/Data
							src/Data/NiiLoader
							src/Alg
							src/Alg/ICP
							src/Alg/NonRigid
							src/Alg/Morphing
							Extern/
							Extern/eigen
							)
FILE( GLOB Project_SRCS "src/*.*" 
						"src/*/*.*"
						"src/*/*/*.*"
						)

# SET QT UI VARS
SET( MainWin_UIS 	src/UI/CranioViewer.ui
					src/UI/MorphingViewer.ui 
					src/UI/DistMapConfigDialog.ui
					)
SET( MainWin_MOC_HDRS 	src/UI/CranioViewer.h
						src/UI/MorphingViewer.h
						src/UI/ImageSliceWidget.h
						src/UI/TrWidget.h 
						src/UI/MainViewer.h 
						src/UI/VTKEventHandler.h
						src/UI/DistMapConfigDialog.h
						src/Alg/ICP/ICPWrapper.h 
						src/Alg/NonRigid/NonRigidWrapper.h
						)

# QT UI WRAPPER
QT5_WRAP_UI( MainWin_UIS_H 		${MainWin_UIS} )
QT5_WRAP_CPP( MainWin_MOC_SRCS 	${MainWin_MOC_HDRS} )

# SET PREDEFINES
ADD_DEFINITIONS( -DNOMINMAX )
ADD_DEFINITIONS( -D_USE_MATH_DEFINES )
ADD_DEFINITIONS( -D_CRT_SECURE_NO_WARNINGS )


#SET( QGLViewer_lib_debug ${CMAKE_CURRENT_SOURCE_DIR}/lib/QGLViewerd2.lib )
#SET( QGLViewer_lib_release ${CMAKE_CURRENT_SOURCE_DIR}/lib/QGLViewer2.lib )

# SET EXECUTABLE RELATED
INCLUDE_DIRECTORIES( ${Project_INCLUDE_DIR}
					 ${CMAKE_CURRENT_BINARY_DIR}
					 ${BOOST_DIR}
					 )

ADD_EXECUTABLE( ${PROJECT_NAME} 
				${Project_SRCS}  
				${MainWin_UIS_H} 
				${MainWin_MOC_SRCS}
				${SVD_SRCS} 
				)
				
TARGET_LINK_LIBRARIES( ${PROJECT_NAME} 
					   ${Qt5Widgets_LIBRARIES} 
					   ${Qt5Xml_LIBRARIES}
					   ${Qt5Core_LIBRARIES}
					   #${Qt5OpenGL_LIBRARIES}
					   #${OPENGL_LIBRARIES}
					   #${OpenCV_LIBS}
					   #${GLEW_LIBRARIES}
					   #${ITK_LIBRARIES}
					   ${vtk_lib_renew}
					   #debug ${QGLViewer_lib_debug} 
					   #optimized ${QGLViewer_lib_release} 
					   )	

# SET PROJECT FILE GROUPS
FILE( GLOB UI 				"src/UI/*.*" )
FILE( GLOB Data_NiiLoader 	"src/Data/NiiLoader/*.*" )
FILE( GLOB Data				"src/Data/*.*" )
FILE( GLOB Alg_ICP 			"src/Alg/ICP/*.*" )
FILE( GLOB Alg_NonRigid		"src/Alg/NonRigid/*.*" )
FILE( GLOB Alg_Morphing		"src/Alg/Morphing/*.*" )
					   
SOURCE_GROUP( "UI" 				FILES ${UI} )
SOURCE_GROUP( "Data" 			FILES ${Data} )
SOURCE_GROUP( "Data\\Nii" 		FILES ${Data_NiiLoader} )
SOURCE_GROUP( "Alg\\ICP" 		FILES ${Alg_ICP} )
SOURCE_GROUP( "Alg\\NonRigid" 	FILES ${Alg_NonRigid} )
SOURCE_GROUP( "Alg\\Morphing" 	FILES ${Alg_Morphing} )

#MESSAGE( STATUS "${ITK_LIBRARIES}" )
MESSAGE( STATUS "${VTK_LIBRARIES}" )
MESSAGE( STATUS "${vtk_lib_renew}" )
MESSAGE( STATUS "${VTK_INCLUDE_DIRECTORIES}" )
