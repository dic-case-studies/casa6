//# QtDBusXmlApp.qo.h: Abstract parent to use the CASA DBus server.
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
#ifndef QTDBUSXMLAPP_QO_H_
#define QTDBUSXMLAPP_QO_H_

#include <casaqt/QtUtilities/QtDBusApp.h>
#include <casaqt/QtUtilities/QtDBusXML.h>

#include <QDBusAbstractAdaptor>
#include <QDBusConnection>

namespace casa {

//# Forward Declarations.
class QtDBusXmlAppAdaptor;


// Abstract parent of classes that want to register themselves with CASA's
// DBus server.  The class also provides public static methods for other code
// to send DBus messages without needing to subclass and register.  This class
// hides the details of how the DBus communication happens behind the scenes
// and presents an interface using CASA objects like Records.  Applications
// that want to use this DBus communication layer need to know about the
// following for each method call:
// * The object name.  Each object that registers with CASA's DBus server must
//   do so with a unique name.  Outside objects then can send messages to that
//   object by using its name.  Names have some restrictions; see
//   dbusRegisterSelf().  You can check if a certain name is currently
//   registered using dbusNameIsRegistered().
// * The method name.  Each object keeps track of what methods it supports, and
//   this process happens dynamically on a per-message basis.
// * Whether to call asynchronously or not.  See dbusXmlCallNoRet().
// * The method parameters.  The parameters are set using name/value pairs in a
//   Record.  Like the method name, this is on a dynamic, per-object basis.
//   IMPORTANT: not all types are currently supported; see QtDBusXML
//   documentation to see what is currently implemented.
// * The method return value.  IMPORTANT: not all types are currently
//   supported; see QtDBusXML documentation to see what is currently
//   implemented.
//
// If a class wants to receive DBus communications, it must take the following
// steps:
// 1) Subclass QtDBusXmlApp.
// 2) Implement the dbusRunXmlMethod() method, which is called when a DBus method
//    call is received.  For a discussion of the parameters to this method, see
//    above.
// 3) Call dbusRegisterSelf() with a unique name.
//
// Classes that wish to send but not receive messages do NOT need to subclass
// QtDBusXmlApp and can just use the public static methods to call methods on
// DBus-registered objects.  For a discussion of the parameters to these
// static methods, see above.
class QtDBusXmlApp : public QtDBusApp {
    
    //# Friend class declarations.
    friend class QtDBusXmlAppAdaptor;

public:
    
    // Calls the given method on the object with the given name that is
    // registered with CASA's DBus server, using the given parameters.  The
    // given from name is sent to the remote object, but not used otherwise. If
    // isAsync is true, then the remote method runs asynchronously, which means
    // that control returns immediately after sending the message.  This method
    // does NOT give a return value, even if the remote method does.  Returns
    // true for success, false for failure.  Will fail if the given object
    // name is not registered with CASA's DBus server.
    static bool dbusXmlCallNoRet(const casacore::String& fromName,
            const casacore::String& objectName, const casacore::String& methodName,
            const casacore::Record& parameters, bool isAsync = false);
    
    // Like dbusXmlCallNoRet(), except that if the remote method has a
    // return value of the given type, then the value is set accordingly.  If
    // there is no return value or it is a different type, the value is not
    // set.  If retValueSet is given, it will be set to true if the return
    // value was set and false otherwise.
    // <group>
    static bool dbusXmlCall(const casacore::String& fromName,
            const casacore::String& objectName, const casacore::String& methodName,
            const casacore::Record& parameters, bool& retValue,
            bool* retValueSet = NULL);
    static bool dbusXmlCall(const casacore::String& fromName,
            const casacore::String& objectName, const casacore::String& methodName,
            const casacore::Record& parameters, int& retValue,
            bool* retValueSet = NULL);
    static bool dbusXmlCall(const casacore::String& fromName,
            const casacore::String& objectName, const casacore::String& methodName,
            const casacore::Record& parameters, casacore::uInt& retValue,
            bool* retValueSet = NULL);
    static bool dbusXmlCall(const casacore::String& fromName,
            const casacore::String& objectName, const casacore::String& methodName,
            const casacore::Record& parameters, double& retValue,
            bool* retValueSet = NULL);
    static bool dbusXmlCall(const casacore::String& fromName,
            const casacore::String& objectName, const casacore::String& methodName,
            const casacore::Record& parameters, casacore::String& retValue,
            bool* retValueSet = NULL);
    static bool dbusXmlCall(const casacore::String& fromName,
            const casacore::String& objectName, const casacore::String& methodName,
            const casacore::Record& parameters, casacore::Record& retValue,
            bool* retValueSet = NULL);
    static bool dbusXmlCall(const casacore::String& fromName,
            const casacore::String& objectName, const casacore::String& methodName,
            const casacore::Record& parameters, casacore::Array<bool>& retValue,
            bool* retValueSet = NULL);
    static bool dbusXmlCall(const casacore::String& fromName,
            const casacore::String& objectName, const casacore::String& methodName,
            const casacore::Record& parameters, casacore::Array<int>& retValue,
            bool* retValueSet = NULL);
    // </group>
    
protected:
    // Constructor.
    QtDBusXmlApp();
    
    // Destructor.  Unregisters from the CASA DBus server if needed.
    virtual ~QtDBusXmlApp();
    
    
    // ABSTRACT METHODS //
    
    // Runs the method with the specified name using the given parameters and
    // placing the return value, if any, into the given retValue record.  NOTE:
    // when defining the return value, the name doesn't matter because the
    // first entry is used.  The caller name, and whether this is an
    // asynchronous call or not, are also provided but do not need to be used.
    // Note, however, that asynchronous method calls will NOT use a return
    // value even if one is set.
    virtual void dbusRunXmlMethod(const casacore::String& methodName,
            const casacore::Record& parameters, casacore::Record& retValue,
            const casacore::String& callerName, bool isAsync) = 0;
    
    
    // VIRTUAL METHODS //
    
    // Method that can be overridden if the subclass wants to be informed
    // whenever ANY dbus message is received, even if this object is not the
    // intended recipient.  Note that most applications won't need to do this
    // (and probably shouldn't) since dbusRunXmlMethod() will be called with the
    // appropriate parameters if this object is the intended recipient.
    virtual void dbusXmlReceived(const QtDBusXML& xml) { (void)xml; }
    
    
    // IMPLEMENTED METHODS //
    
    // Registers this object with CASA's DBus server, if it is not already,
    // with the given name and returns whether or not the registration
    // succeeded.  If the name is blank, then the last set name is used, UNLESS
    // this is the first time registering in which case the registration will
    // fail.  The registration name MUST contain only numbers, letters, and
    // underscores, and MUST be unique for the DBus server; trying to register
    // with a name that is already in use will result in registration failure.
    // (Note, however, that the name needs ONLY to be unique within the CASA
    // DBus application names rather than all system-wide DBus application
    // names.)  Trying to register when already registered (see
    // dbusSelfIsRegistered()) will result in registration failure; to change
    // names, you must unregister and then reregister with the new name.
    bool dbusRegisterSelf(const casacore::String& name = "");
    
    // Unregisters this object with CASA's DBus server, if it is registered.
    void dbusUnregisterSelf();
    
    // Returns true if this object is currently registered with CASA's DBus
    // server, false otherwise.
    bool dbusSelfIsRegistered() const;
    
    // Returns the name that this object is registered with with CASA's DBus
    // server, or an empty casacore::String if this application is not currently
    // registered.
    casacore::String dbusSelfRegisteredName() const;
    
    // Calls the static version of the method with this application's name.
    // <group>
    bool dbusXmlCallNoRet(const casacore::String& objectName,
            const casacore::String& methodName, const casacore::Record& parameters,
            bool isAsync = false) {
        return dbusXmlCallNoRet(dbusSelfRegisteredName(), objectName,
                                   methodName, parameters, isAsync); }
    bool dbusXmlCall(const casacore::String& objectName, const casacore::String& methodName,
            const casacore::Record& parameters, bool& retValue,
            bool* retValueSet = NULL) {
        return dbusXmlCall(dbusSelfRegisteredName(), objectName, methodName,
                              parameters, retValue, retValueSet); }
    bool dbusXmlCall(const casacore::String& objectName, const casacore::String& methodName,
            const casacore::Record& parameters, int& retValue,
            bool* retValueSet = NULL) {
        return dbusXmlCall(dbusSelfRegisteredName(), objectName, methodName,
                              parameters, retValue, retValueSet); }
    bool dbusXmlCall(const casacore::String& objectName, const casacore::String& methodName,
            const casacore::Record& parameters, casacore::uInt& retValue,
            bool* retValueSet = NULL) {
        return dbusXmlCall(dbusSelfRegisteredName(), objectName, methodName,
                              parameters, retValue, retValueSet); }
    bool dbusXmlCall(const casacore::String& objectName, const casacore::String& methodName,
            const casacore::Record& parameters, double& retValue,
            bool* retValueSet = NULL) {
        return dbusXmlCall(dbusSelfRegisteredName(), objectName, methodName,
                              parameters, retValue, retValueSet); }
    bool dbusXmlCall(const casacore::String& objectName, const casacore::String& methodName,
            const casacore::Record& parameters, casacore::String& retValue,
            bool* retValueSet = NULL) {
        return dbusXmlCall(dbusSelfRegisteredName(), objectName, methodName,
                              parameters, retValue, retValueSet); }
    bool dbusXmlCall(const casacore::String& objectName, const casacore::String& methodName,
            const casacore::Record& parameters, casacore::Record& retValue,
            bool* retValueSet = NULL) {
        return dbusXmlCall(dbusSelfRegisteredName(), objectName, methodName,
                              parameters, retValue, retValueSet); }
    bool dbusXmlCall(const casacore::String& objectName, const casacore::String& methodName,
            const casacore::Record& parameters, casacore::Array<bool>& retValue,
            bool* retValueSet = NULL) {
        return dbusXmlCall(dbusSelfRegisteredName(), objectName, methodName,
                              parameters, retValue, retValueSet); }
    bool dbusXmlCall(const casacore::String& objectName, const casacore::String& methodName,
            const casacore::Record& parameters, casacore::Array<int>& retValue,
            bool* retValueSet = NULL) {
        return dbusXmlCall(dbusSelfRegisteredName(), objectName, methodName,
                              parameters, retValue, retValueSet); }
    // </group>
    
private:
    
    // Flag for whether the application is currently registered or not.
    bool dbusRegistered_;
    
    // Name that the application is registered with.
    QString dbusName_;
    
    // DBus adaptor.
    QtDBusXmlAppAdaptor* dbusAdaptor_;

    
    // Method for when one of the slots in the adaptor is activated.  First
    // sends to dbusXmlReceived(), then to dbusRunXmlMethod().
    void dbusSlot(QtDBusXML& xml);
    
    
    // Private Static Methods //
    
    // Helper method for calling remote methods.
    static bool dbusXmlCall(const casacore::String& from, const casacore::String& to,
            const casacore::String& methodName, bool methodIsAsync,
            const casacore::Record& parameters, casacore::Record* retValue);
};


// Subclass of QDBusAbstractAdaptor for use with CASA's QtDBusXmlApp class.  This
// class is a very thin layer on top of QtDBusXmlApp.
class QtDBusXmlAppAdaptor : public QDBusAbstractAdaptor {
     Q_OBJECT
    
    // Interface name definition.
    // <group>
    #define CASA_DBUS_XML_INTERFACE "edu.nrao.casa.QtDBusXmlApp"
    Q_CLASSINFO("D-Bus Interface", "edu.nrao.casa.QtDBusXmlApp")
    // </group>
    
    //# Friend class declarations.
    friend class QtDBusXmlApp;
    
public slots:
    // Slot for receiving messages.  If its name is changed,
    // QtDBusXmlApp::DBUS_MESSAGE_SLOT must be updated.
    QString xmlSlot(const QString& xml);
    
private:
    // Constructor which takes the application.
    QtDBusXmlAppAdaptor(QtDBusXmlApp& app);
    
    // Destructor.
    ~QtDBusXmlAppAdaptor();
    
    
    // Application.
    QtDBusXmlApp& itsApp_;
    
    // Dummy object.
    QObject* itsObject_;
};

}

#endif /* QTDBUSAPP_QO_H_ */
