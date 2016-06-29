/****************************************************************************
** Meta object code from reading C++ file 'MainWindow.h'
**
** Created: Thu Feb 28 13:51:58 2013
**      by: The Qt Meta Object Compiler version 62 (Qt 4.6.3)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "MainWindow.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'MainWindow.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 62
#error "This file was generated using the moc from 4.6.3. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_MainWindow[] = {

 // content:
       4,       // revision
       0,       // classname
       0,    0, // classinfo
      18,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: signature, parameters, type, tag, flags
      12,   11,   11,   11, 0x08,
      22,   11,   11,   11, 0x08,
      30,   11,   11,   11, 0x08,
      37,   11,   11,   11, 0x08,
      52,   11,   11,   11, 0x08,
      66,   11,   11,   11, 0x08,
      79,   11,   11,   11, 0x08,
      87,   11,   11,   11, 0x08,
     106,   96,   11,   11, 0x08,
     127,   96,   11,   11, 0x08,
     148,   96,   11,   11, 0x08,
     160,   96,   11,   11, 0x08,
     172,   96,   11,   11, 0x08,
     192,   96,   11,   11, 0x08,
     208,   96,   11,   11, 0x08,
     224,   96,   11,   11, 0x08,
     242,   96,   11,   11, 0x08,
     261,   96,   11,   11, 0x08,

       0        // eod
};

static const char qt_meta_stringdata_MainWindow[] = {
    "MainWindow\0\0iterate()\0start()\0stop()\0"
    "continue_sim()\0start_video()\0stop_video()\0"
    "reset()\0finish()\0rate0_int\0"
    "set_reprodrate0(int)\0set_reprodrate1(int)\0"
    "set_r0(int)\0set_r1(int)\0set_dediffrate(int)\0"
    "set_delta0(int)\0set_delta1(int)\0"
    "set_shedrate(int)\0set_stratrate(int)\0"
    "set_directstratfrac(int)\0"
};

const QMetaObject MainWindow::staticMetaObject = {
    { &QDialog::staticMetaObject, qt_meta_stringdata_MainWindow,
      qt_meta_data_MainWindow, 0 }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &MainWindow::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *MainWindow::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *MainWindow::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_MainWindow))
        return static_cast<void*>(const_cast< MainWindow*>(this));
    return QDialog::qt_metacast(_clname);
}

int MainWindow::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QDialog::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: iterate(); break;
        case 1: start(); break;
        case 2: stop(); break;
        case 3: continue_sim(); break;
        case 4: start_video(); break;
        case 5: stop_video(); break;
        case 6: reset(); break;
        case 7: finish(); break;
        case 8: set_reprodrate0((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 9: set_reprodrate1((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 10: set_r0((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 11: set_r1((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 12: set_dediffrate((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 13: set_delta0((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 14: set_delta1((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 15: set_shedrate((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 16: set_stratrate((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 17: set_directstratfrac((*reinterpret_cast< int(*)>(_a[1]))); break;
        default: ;
        }
        _id -= 18;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
