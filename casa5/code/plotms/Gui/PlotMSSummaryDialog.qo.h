//# Copyright (C) 2009
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
//# $Id: $
#ifndef PLOTMSSUMMARYDIALOG_QO_H
#define PLOTMSSUMMARYDIALOG_QO_H

#include <QDialog>
#include <ui/ui_PlotMSSummaryDialog.h>
#include <plotms/PlotMS/PlotMSConstants.h>
#include <casa/BasicSL/String.h>

namespace casa {

/**
 * Contains controls for customizing and producing a summary of the
 * plot.
 */


class PlotMSSummaryDialog : public QDialog {
    Q_OBJECT

public:
    PlotMSSummaryDialog(QDialog *parent = 0);

    //Is the summary verbose?
    bool isVerbose() const;
    casacore::String getFileName() const;
    void filesChanged(const std::vector<casacore::String>& fileNamees);
    //Return the summary type.
    PMS::SummaryType getSummaryType() const;
    PMS::CTSummaryType getCTSummaryType() const;
	inline bool isMS() { return isMS_; }
    ~PlotMSSummaryDialog();

private slots:
	void closeDialog();
	void summarize();

private:
    PlotMSSummaryDialog* summarizeDialog;
    Ui::PlotMSSummaryDialogClass ui;

	bool isMS_;
};

}

#endif // PLOTMSSUMMARYDIALOG_QO_H
