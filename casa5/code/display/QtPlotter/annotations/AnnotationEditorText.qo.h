#ifndef ANNOTATIONEDITORTEXT_QO_H
#define ANNOTATIONEDITORTEXT_QO_H

#include <QWidget>
#include <ui/ui_AnnotationEditorText.h>

namespace casa {
	class AnnotationEditorText : public QWidget {
		Q_OBJECT

	public:
		AnnotationEditorText(QWidget *parent = 0);
		QString getLabel() const;
		QFont getFontFamily() const;
		bool isBold() const;
		bool isItalic() const;
		int getFontSize() const;
		~AnnotationEditorText();

	private:
		Ui::AnnotationEditorTextClass ui;
	};
}
#endif // ANNOTATIONEDITORTEXT_QO_H
