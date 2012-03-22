/********************************************************************************
** Form generated from reading UI file 'finestra.ui'
**
** Created: Mon Feb 6 00:06:49 2012
**      by: Qt User Interface Compiler version 4.7.4
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_FINESTRA_H
#define UI_FINESTRA_H

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QCheckBox>
#include <QtGui/QGroupBox>
#include <QtGui/QHBoxLayout>
#include <QtGui/QHeaderView>
#include <QtGui/QLabel>
#include <QtGui/QMainWindow>
#include <QtGui/QMenuBar>
#include <QtGui/QProgressBar>
#include <QtGui/QPushButton>
#include <QtGui/QRadioButton>
#include <QtGui/QSpacerItem>
#include <QtGui/QSpinBox>
#include <QtGui/QStatusBar>
#include <QtGui/QVBoxLayout>
#include <QtGui/QWidget>

QT_BEGIN_NAMESPACE

class Ui_MainWindow
{
public:
    QWidget *centralwidget;
    QWidget *widget;
    QHBoxLayout *horizontalLayout_5;
    QProgressBar *progressBar;
    QVBoxLayout *verticalLayout_4;
    QPushButton *startButton;
    QPushButton *closeButton;
    QWidget *widget1;
    QHBoxLayout *horizontalLayout;
    QGroupBox *groupBox;
    QRadioButton *fromRandom;
    QRadioButton *fromFile;
    QRadioButton *fromIsing;
    QGroupBox *groupBox_2;
    QWidget *widget2;
    QHBoxLayout *horizontalLayout_2;
    QLabel *label;
    QSpinBox *spinBox;
    QWidget *widget3;
    QHBoxLayout *horizontalLayout_3;
    QLabel *label_2;
    QSpinBox *spinBox_3;
    QWidget *widget4;
    QHBoxLayout *horizontalLayout_4;
    QLabel *label_3;
    QSpinBox *spinBox_2;
    QWidget *widget5;
    QVBoxLayout *verticalLayout_2;
    QSpacerItem *verticalSpacer;
    QVBoxLayout *verticalLayout;
    QCheckBox *checkBox;
    QCheckBox *checkBox_2;
    QCheckBox *checkBox_3;
    QCheckBox *checkBox_4;
    QGroupBox *groupBox_3;
    QWidget *widget6;
    QVBoxLayout *verticalLayout_3;
    QCheckBox *checkRohlin;
    QCheckBox *checkRohlinTop;
    QCheckBox *checkRidotta;
    QCheckBox *checkRidottaTop;
    QCheckBox *checkFuzzy;
    QCheckBox *checkBox_5;
    QCheckBox *checkSomiglianza;
    QMenuBar *menubar;
    QStatusBar *statusbar;

    void setupUi(QMainWindow *MainWindow)
    {
        if (MainWindow->objectName().isEmpty())
            MainWindow->setObjectName(QString::fromUtf8("MainWindow"));
        MainWindow->resize(614, 448);
        QSizePolicy sizePolicy(QSizePolicy::Minimum, QSizePolicy::Minimum);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(MainWindow->sizePolicy().hasHeightForWidth());
        MainWindow->setSizePolicy(sizePolicy);
        centralwidget = new QWidget(MainWindow);
        centralwidget->setObjectName(QString::fromUtf8("centralwidget"));
        widget = new QWidget(centralwidget);
        widget->setObjectName(QString::fromUtf8("widget"));
        widget->setGeometry(QRect(10, 340, 601, 62));
        horizontalLayout_5 = new QHBoxLayout(widget);
        horizontalLayout_5->setObjectName(QString::fromUtf8("horizontalLayout_5"));
        horizontalLayout_5->setContentsMargins(0, 0, 0, 0);
        progressBar = new QProgressBar(widget);
        progressBar->setObjectName(QString::fromUtf8("progressBar"));
        progressBar->setEnabled(false);
        progressBar->setValue(0);

        horizontalLayout_5->addWidget(progressBar);

        verticalLayout_4 = new QVBoxLayout();
        verticalLayout_4->setObjectName(QString::fromUtf8("verticalLayout_4"));
        startButton = new QPushButton(widget);
        startButton->setObjectName(QString::fromUtf8("startButton"));

        verticalLayout_4->addWidget(startButton);

        closeButton = new QPushButton(widget);
        closeButton->setObjectName(QString::fromUtf8("closeButton"));

        verticalLayout_4->addWidget(closeButton);


        horizontalLayout_5->addLayout(verticalLayout_4);

        widget1 = new QWidget(centralwidget);
        widget1->setObjectName(QString::fromUtf8("widget1"));
        widget1->setGeometry(QRect(0, 10, 611, 331));
        horizontalLayout = new QHBoxLayout(widget1);
        horizontalLayout->setObjectName(QString::fromUtf8("horizontalLayout"));
        horizontalLayout->setContentsMargins(0, 0, 0, 0);
        groupBox = new QGroupBox(widget1);
        groupBox->setObjectName(QString::fromUtf8("groupBox"));
        fromRandom = new QRadioButton(groupBox);
        fromRandom->setObjectName(QString::fromUtf8("fromRandom"));
        fromRandom->setGeometry(QRect(10, 20, 97, 20));
        fromFile = new QRadioButton(groupBox);
        fromFile->setObjectName(QString::fromUtf8("fromFile"));
        fromFile->setGeometry(QRect(10, 140, 97, 20));
        fromFile->setChecked(true);
        fromIsing = new QRadioButton(groupBox);
        fromIsing->setObjectName(QString::fromUtf8("fromIsing"));
        fromIsing->setGeometry(QRect(10, 250, 97, 20));
        fromRandom->raise();
        fromFile->raise();
        fromIsing->raise();

        horizontalLayout->addWidget(groupBox);

        groupBox_2 = new QGroupBox(widget1);
        groupBox_2->setObjectName(QString::fromUtf8("groupBox_2"));
        widget2 = new QWidget(groupBox_2);
        widget2->setObjectName(QString::fromUtf8("widget2"));
        widget2->setGeometry(QRect(10, 30, 191, 45));
        horizontalLayout_2 = new QHBoxLayout(widget2);
        horizontalLayout_2->setObjectName(QString::fromUtf8("horizontalLayout_2"));
        horizontalLayout_2->setContentsMargins(0, 0, 0, 0);
        label = new QLabel(widget2);
        label->setObjectName(QString::fromUtf8("label"));
        label->setWordWrap(true);

        horizontalLayout_2->addWidget(label);

        spinBox = new QSpinBox(widget2);
        spinBox->setObjectName(QString::fromUtf8("spinBox"));
        spinBox->setMinimum(1);
        spinBox->setMaximum(5000);
        spinBox->setValue(1000);

        horizontalLayout_2->addWidget(spinBox);

        widget3 = new QWidget(groupBox_2);
        widget3->setObjectName(QString::fromUtf8("widget3"));
        widget3->setGeometry(QRect(10, 80, 191, 41));
        horizontalLayout_3 = new QHBoxLayout(widget3);
        horizontalLayout_3->setObjectName(QString::fromUtf8("horizontalLayout_3"));
        horizontalLayout_3->setContentsMargins(0, 0, 0, 0);
        label_2 = new QLabel(widget3);
        label_2->setObjectName(QString::fromUtf8("label_2"));
        label_2->setWordWrap(true);

        horizontalLayout_3->addWidget(label_2);

        spinBox_3 = new QSpinBox(widget3);
        spinBox_3->setObjectName(QString::fromUtf8("spinBox_3"));
        QSizePolicy sizePolicy1(QSizePolicy::Fixed, QSizePolicy::Fixed);
        sizePolicy1.setHorizontalStretch(0);
        sizePolicy1.setVerticalStretch(0);
        sizePolicy1.setHeightForWidth(spinBox_3->sizePolicy().hasHeightForWidth());
        spinBox_3->setSizePolicy(sizePolicy1);
        spinBox_3->setMinimum(1);
        spinBox_3->setMaximum(5000);
        spinBox_3->setValue(1000);

        horizontalLayout_3->addWidget(spinBox_3);

        widget4 = new QWidget(groupBox_2);
        widget4->setObjectName(QString::fromUtf8("widget4"));
        widget4->setGeometry(QRect(10, 120, 191, 27));
        horizontalLayout_4 = new QHBoxLayout(widget4);
        horizontalLayout_4->setObjectName(QString::fromUtf8("horizontalLayout_4"));
        horizontalLayout_4->setContentsMargins(0, 0, 0, 0);
        label_3 = new QLabel(widget4);
        label_3->setObjectName(QString::fromUtf8("label_3"));
        label_3->setWordWrap(true);

        horizontalLayout_4->addWidget(label_3);

        spinBox_2 = new QSpinBox(widget4);
        spinBox_2->setObjectName(QString::fromUtf8("spinBox_2"));
        spinBox_2->setMaximum(100);
        spinBox_2->setValue(5);

        horizontalLayout_4->addWidget(spinBox_2);

        widget5 = new QWidget(groupBox_2);
        widget5->setObjectName(QString::fromUtf8("widget5"));
        widget5->setGeometry(QRect(10, 150, 164, 152));
        verticalLayout_2 = new QVBoxLayout(widget5);
        verticalLayout_2->setObjectName(QString::fromUtf8("verticalLayout_2"));
        verticalLayout_2->setContentsMargins(0, 0, 0, 0);
        verticalSpacer = new QSpacerItem(20, 40, QSizePolicy::Minimum, QSizePolicy::Expanding);

        verticalLayout_2->addItem(verticalSpacer);

        verticalLayout = new QVBoxLayout();
        verticalLayout->setObjectName(QString::fromUtf8("verticalLayout"));
        checkBox = new QCheckBox(widget5);
        checkBox->setObjectName(QString::fromUtf8("checkBox"));

        verticalLayout->addWidget(checkBox);

        checkBox_2 = new QCheckBox(widget5);
        checkBox_2->setObjectName(QString::fromUtf8("checkBox_2"));

        verticalLayout->addWidget(checkBox_2);

        checkBox_3 = new QCheckBox(widget5);
        checkBox_3->setObjectName(QString::fromUtf8("checkBox_3"));

        verticalLayout->addWidget(checkBox_3);

        checkBox_4 = new QCheckBox(widget5);
        checkBox_4->setObjectName(QString::fromUtf8("checkBox_4"));

        verticalLayout->addWidget(checkBox_4);


        verticalLayout_2->addLayout(verticalLayout);


        horizontalLayout->addWidget(groupBox_2);

        groupBox_3 = new QGroupBox(widget1);
        groupBox_3->setObjectName(QString::fromUtf8("groupBox_3"));
        widget6 = new QWidget(groupBox_3);
        widget6->setObjectName(QString::fromUtf8("widget6"));
        widget6->setGeometry(QRect(20, 30, 161, 271));
        verticalLayout_3 = new QVBoxLayout(widget6);
        verticalLayout_3->setObjectName(QString::fromUtf8("verticalLayout_3"));
        verticalLayout_3->setContentsMargins(0, 0, 0, 0);
        checkRohlin = new QCheckBox(widget6);
        checkRohlin->setObjectName(QString::fromUtf8("checkRohlin"));

        verticalLayout_3->addWidget(checkRohlin);

        checkRohlinTop = new QCheckBox(widget6);
        checkRohlinTop->setObjectName(QString::fromUtf8("checkRohlinTop"));
        checkRohlinTop->setEnabled(false);

        verticalLayout_3->addWidget(checkRohlinTop);

        checkRidotta = new QCheckBox(widget6);
        checkRidotta->setObjectName(QString::fromUtf8("checkRidotta"));
        checkRidotta->setEnabled(false);

        verticalLayout_3->addWidget(checkRidotta);

        checkRidottaTop = new QCheckBox(widget6);
        checkRidottaTop->setObjectName(QString::fromUtf8("checkRidottaTop"));
        checkRidottaTop->setEnabled(false);

        verticalLayout_3->addWidget(checkRidottaTop);

        checkFuzzy = new QCheckBox(widget6);
        checkFuzzy->setObjectName(QString::fromUtf8("checkFuzzy"));

        verticalLayout_3->addWidget(checkFuzzy);

        checkBox_5 = new QCheckBox(widget6);
        checkBox_5->setObjectName(QString::fromUtf8("checkBox_5"));
        checkBox_5->setEnabled(false);

        verticalLayout_3->addWidget(checkBox_5);

        checkSomiglianza = new QCheckBox(widget6);
        checkSomiglianza->setObjectName(QString::fromUtf8("checkSomiglianza"));

        verticalLayout_3->addWidget(checkSomiglianza);


        horizontalLayout->addWidget(groupBox_3);

        MainWindow->setCentralWidget(centralwidget);
        startButton->raise();
        groupBox->raise();
        groupBox_2->raise();
        groupBox_3->raise();
        progressBar->raise();
        closeButton->raise();
        menubar = new QMenuBar(MainWindow);
        menubar->setObjectName(QString::fromUtf8("menubar"));
        menubar->setGeometry(QRect(0, 0, 614, 23));
        MainWindow->setMenuBar(menubar);
        statusbar = new QStatusBar(MainWindow);
        statusbar->setObjectName(QString::fromUtf8("statusbar"));
        MainWindow->setStatusBar(statusbar);
#ifndef QT_NO_SHORTCUT
        label->setBuddy(spinBox);
        label_2->setBuddy(spinBox_3);
        label_3->setBuddy(spinBox_2);
#endif // QT_NO_SHORTCUT

        retranslateUi(MainWindow);
        QObject::connect(closeButton, SIGNAL(clicked()), MainWindow, SLOT(close()));
        QObject::connect(checkBox, SIGNAL(clicked(bool)), checkRohlin, SLOT(toggle()));
        QObject::connect(checkBox, SIGNAL(clicked(bool)), checkRidotta, SLOT(setChecked(bool)));
        QObject::connect(checkBox, SIGNAL(clicked(bool)), checkFuzzy, SLOT(setChecked(bool)));
        QObject::connect(checkBox, SIGNAL(clicked(bool)), checkSomiglianza, SLOT(setChecked(bool)));
        QObject::connect(checkRidotta, SIGNAL(clicked(bool)), checkRidottaTop, SLOT(setChecked(bool)));
        QObject::connect(checkRohlin, SIGNAL(clicked(bool)), checkRohlinTop, SLOT(setChecked(bool)));
        QObject::connect(checkRohlin, SIGNAL(clicked(bool)), checkRidotta, SLOT(setEnabled(bool)));
        QObject::connect(checkFuzzy, SIGNAL(clicked(bool)), checkBox_5, SLOT(setEnabled(bool)));
        QObject::connect(checkRohlin, SIGNAL(clicked(bool)), checkRohlinTop, SLOT(setEnabled(bool)));

        QMetaObject::connectSlotsByName(MainWindow);
    } // setupUi

    void retranslateUi(QMainWindow *MainWindow)
    {
        MainWindow->setWindowTitle(QApplication::translate("MainWindow", "Calcola distanze tra sequenze", 0, QApplication::UnicodeUTF8));
        startButton->setText(QApplication::translate("MainWindow", "Esegui", 0, QApplication::UnicodeUTF8));
        closeButton->setText(QApplication::translate("MainWindow", "Chiudi", 0, QApplication::UnicodeUTF8));
        groupBox->setTitle(QApplication::translate("MainWindow", "Sorgente", 0, QApplication::UnicodeUTF8));
        fromRandom->setText(QApplication::translate("MainWindow", "Random", 0, QApplication::UnicodeUTF8));
        fromFile->setText(QApplication::translate("MainWindow", "File", 0, QApplication::UnicodeUTF8));
        fromIsing->setText(QApplication::translate("MainWindow", "Ising", 0, QApplication::UnicodeUTF8));
        groupBox_2->setTitle(QApplication::translate("MainWindow", "Parametri", 0, QApplication::UnicodeUTF8));
        label->setText(QApplication::translate("MainWindow", "Numero di sequenze", 0, QApplication::UnicodeUTF8));
        label_2->setText(QApplication::translate("MainWindow", "Lunghezza massima sequenze", 0, QApplication::UnicodeUTF8));
        label_3->setText(QApplication::translate("MainWindow", "Salto massimo", 0, QApplication::UnicodeUTF8));
        checkBox->setText(QApplication::translate("MainWindow", "Non calcolare distanze", 0, QApplication::UnicodeUTF8));
        checkBox_2->setText(QApplication::translate("MainWindow", "Verboso", 0, QApplication::UnicodeUTF8));
        checkBox_3->setText(QApplication::translate("MainWindow", "Scrivi risultati", 0, QApplication::UnicodeUTF8));
        checkBox_4->setText(QApplication::translate("MainWindow", "Riduci alfabeto", 0, QApplication::UnicodeUTF8));
        groupBox_3->setTitle(QApplication::translate("MainWindow", "Distanze da calcolare", 0, QApplication::UnicodeUTF8));
        checkRohlin->setText(QApplication::translate("MainWindow", "Rohlin", 0, QApplication::UnicodeUTF8));
        checkRohlinTop->setText(QApplication::translate("MainWindow", "Topologica", 0, QApplication::UnicodeUTF8));
        checkRidotta->setText(QApplication::translate("MainWindow", "Rohlin ridotta", 0, QApplication::UnicodeUTF8));
        checkRidottaTop->setText(QApplication::translate("MainWindow", "Topologica", 0, QApplication::UnicodeUTF8));
        checkFuzzy->setText(QApplication::translate("MainWindow", "Fuzzy", 0, QApplication::UnicodeUTF8));
        checkBox_5->setText(QApplication::translate("MainWindow", "Ridotta", 0, QApplication::UnicodeUTF8));
        checkSomiglianza->setText(QApplication::translate("MainWindow", "Somiglianza", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class MainWindow: public Ui_MainWindow {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_FINESTRA_H
