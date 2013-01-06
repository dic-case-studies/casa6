//# Copyright (C) 1994,1995,1996,1997,1998,1999,2000
//# Associated Universities, Inc. Washington DC, USA.
//#
//# This library is free software; you can redistribute it and/or modify it
//# under the terms of the GNU Library General Public License as published by
//# the Free Software Foundation; either version 2 of the License, or (at your
//# option) any later version.
//#
//# This library is distributed in the hope that it will be useful, but WITHOUT
//# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Library General Public
//# License for more details.
//#
//# You should have received a copy of the GNU Library General Public License
//# along with this library; if not, write to the Free Software Foundation,
//# Inc., 675 Massachusetts Ave, Cambridge, MA 02139, USA.
//#
//# Correspondence concerning AIPS++ should be addressed as follows:
//#        Internet email: aips2-request@nrao.edu.
//#        Postal address: AIPS++ Project Office
//#                        National Radio Astronomy Observatory
//#                        520 Edgemont Road
//#                        Charlottesville, VA 22903-2475 USA
//#

#ifndef COLORCOMBODELEGATE_H_
#define COLORCOMBODELEGATE_H_

#include <QItemDelegate>

class QPainter;

namespace casa {

class ColorComboDelegate : public QItemDelegate {
public:
	ColorComboDelegate(QObject* parent = NULL);
	void setSupportedColors( QStringList colorNames );
	QString getNamedColor( int index ) const;
	void paint( QPainter* painter, const QStyleOptionViewItem& option, const QModelIndex& index ) const;
	virtual ~ColorComboDelegate();
private:
	QStringList colorNames;
};

} /* namespace casa */
#endif /* COLORCOMBODELEGATE_H_ */
