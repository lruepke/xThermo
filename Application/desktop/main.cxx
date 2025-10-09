/*
 * Copyright 2007 Sandia Corporation.
 * Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
 * license for use of this work by or on behalf of the
 * U.S. Government. Redistribution and use in source and binary forms, with
 * or without modification, are permitted provided that this Notice and any
 * statement of authorship are reproduced on all copies.
 */

// QT includes 
#include <QApplication>
#include <QSurfaceFormat>

//#include "QVTKOpenGLWidget.h"
#include "QVTKRenderWidget.h"
#include "MainWindow.h"
#include "config.h"
#include "SWEOSbash.h"
// #include "touchbar.h"
#include "vtkAutoInit.h"

//VTK_MODULE_INIT(vtkRenderingContext2D_AUTOINIT)

#include "vtkContextScene.h"
#include "vtkContextView.h"
#include "vtkPiecewiseControlPointsItem.h"
#include "vtkPiecewiseFunction.h"
#include "vtkPiecewiseFunctionItem.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkRenderer.h"
#include "vtkSmartPointer.h"


extern int qInitResources_icons();

int main( int argc, char** argv )
{

  #ifdef _WIN32
      if(argc == 1)
      {
        // needed to ensure appropriate OpenGL context is created for VTK rendering.
            QSurfaceFormat::setDefaultFormat(QVTKOpenGLWidget::defaultFormat());   // must be here
              // QCoreApplication::setAttribute(Qt::AA_EnableHighDpiScaling);
            // QT Stuff
            QApplication app( argc, argv );

            QApplication::setStyle("fusion");

            qInitResources_icons();

            MainWindow myMainWindow;
            //  myMainWindow.showMaximized();
            myMainWindow.show();

            return app.exec();
      }else
      {
        SWEOSbash::bash_run(argc, argv);
      }
  #else
      if(argc==2 || argc==1)
      {
        SWEOSbash::cSWEOSarg arg;
          if(arg.Parse(argc, argv))
          {
            SWEOSbash::bash_run(argc, argv);
          }else
          {
            // needed to ensure appropriate OpenGL context is created for VTK rendering.
            QSurfaceFormat::setDefaultFormat(QVTKRenderWidget::defaultFormat());   // must be here
              // QCoreApplication::setAttribute(Qt::AA_EnableHighDpiScaling);
            // QT Stuff
            QApplication app( argc, argv );

            // {
            //     // Install TouchBarProvider as application delegate
            //     TouchBarProvider *touchBarProvider = [[TouchBarProvider alloc] init];
            //     [touchBarProvider installAsDelegateForApplication:[NSApplication sharedApplication]];
            // }

            QApplication::setStyle("fusion");

            qInitResources_icons();

            MainWindow myMainWindow;
            //  myMainWindow.showMaximized();
            myMainWindow.show();


            // {
            //     // Install TouchBarProvider as window delegate
            //     NSView *view = reinterpret_cast<NSView *>(textEdit.winId());
            //     TouchBarProvider *touchBarProvider = [[TouchBarProvider alloc] init];
            //     [touchBarProvider installAsDelegateForWindow:view.window];
            // }
            return app.exec();
          }
      }else
      {
        SWEOSbash::bash_run(argc, argv);
      }
  #endif
 
 return 0;
}
