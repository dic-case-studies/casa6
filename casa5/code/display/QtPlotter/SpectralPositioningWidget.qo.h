#ifndef SPECTRALPOSITIONINGWIDGET_QO_H
#define SPECTRALPOSITIONINGWIDGET_QO_H

#include <QWidget>
#include <ui/ui_SpectralPositioningWidget.h>
#include <casa/Arrays/Vector.h>
namespace casacore{

	class LogIO;
}

namespace casa {

	class ProfileTaskMonitor;

	class SpectralPositioningWidget : public QWidget {
		Q_OBJECT

	public:
		SpectralPositioningWidget(QWidget *parent = 0);
		void setTaskMonitor( ProfileTaskMonitor* monitor );
		void setLogger( casacore::LogIO* logger );

		void updateRegion( const casacore::Vector<double> px, const casacore::Vector<double> py,
		                   const casacore::Vector<double> wx, const casacore::Vector<double> wy );
		~SpectralPositioningWidget();

	private slots:
		void boxSpecChanged( int index );
		void locationSelectionTypeChanged( int index );
		void locationUnitsChanged( int index );
		void setPosition();

	private:
		void updateUI();
		void updateUIWorldBox();
		void updateUIWorldPoint();
		void updateUIPixelBox();
		void updateUIPixelPoint();
		/**
		  * Initializes the spectrum positioning tab.
		  */
		void initSpectrumPosition();
		void pageUpdate( int selectionIndex, int unitIndex );

		bool populateWorlds( const QList<int> &pixelX, const QList<int> &pixelY,
		                     QList<double> &worldX, QList<double> &worldY );
		bool fillPointWorld( QList<double> &worldX, QList<double> &worldY );
		void fillPointPixel( QList<int> &pixelX, QList<int>&pixelY )const;
		bool fillBoxPixel( QList<int> &pixelX, QList<int>&pixelY );
		bool fillBoxWorld( QList<double> &worldX, QList<double> & worldY );
		bool fillBasedOnBoxSpecification(  const double*  const firstXPix, const double * const firstYPix,
		                                   const double* const secondXPix, const double* const secondYPix,
		                                   double* const blcxPix, double* const blcyPix,
		                                   double* const trcxPix, double* const trcYPix, bool pixels=true );
		double toRadians( bool& valid, QLineEdit * lineEdit );
		void switchBoxLabels( int index, int pageIndex, QLabel* const x1Label, QLabel* const y1Label,
		                      QLabel* const x2Label, QLabel* const y2Label );
		void setPixelLineEdits( double topLeft, double bottomLeft,
		                        double topRight, double bottomRight );
		void setWorldEdits( double topLeft, double bottomLeft,
		                    double topRight, double bottomRight );
		void adjustPoint( const casacore::Vector<double>& newX, const casacore::Vector<double>& newY,
		                  casacore::Vector<double>& xValues, casacore::Vector<double>& yValues );
		Ui::SpectralPositioningWidgetClass ui;

		enum PositionTypeIndex { POINT, BOX, END_POSITION_TYPE };
		enum UnitIndex {RADIAN, PIXEL, END_UNIT };
		QIntValidator* pixelValidator;
		enum StackPages { POINT_PIXEL, POINT_RA_DEC, BOX_PIXEL, BOX_RA_DEC };
		enum BoxSpecificationIndex { TL_LENGTH_HEIGHT, CENTER_LENGTH_HEIGHT, TL_BR, BL_TR,
		                             TL_LENGTH_HEIGHT_WORLD, CENTER_LENGTH_HEIGHT_WORLD, TL_BR_WORLD, BL_TR_WORLD, END_SPEC
		                           };
		QMap<BoxSpecificationIndex,QList<QString> > boxLabelMap;
		ProfileTaskMonitor* profileTaskMonitor;
		casacore::LogIO* logger;
		casacore::Vector<double> pixelXValues;
		casacore::Vector<double> pixelYValues;
		casacore::Vector<double> worldXValues;
		casacore::Vector<double> worldYValues;
	};
}
#endif // SPECTRALPOSITIONINGWIDGET_H
