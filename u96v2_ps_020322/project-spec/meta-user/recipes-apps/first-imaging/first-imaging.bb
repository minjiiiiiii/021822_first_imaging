#
# This file is the first-imaging recipe.
#

BB_DONT_CACHE = "1"
SUMMARY = "Simple first-imaging application"
SECTION = "PETALINUX/apps"
LICENSE = "MIT"
LIC_FILES_CHKSUM = "file://${COMMON_LICENSE_DIR}/MIT;md5=0835ade698e0bcf8506ecda2f7b4f302"

SRC_URI = "file://*.cpp \
		   file://*.h \
		   file://data/* \
           file://Makefile \
		  "

S = "${WORKDIR}"
do_compile[nostamp] = "1"
do_compile() {
		oe_runmake clean
		oe_runmake
}

do_install() {
	     install -d ${D}${bindir}
	     install -m 0755 first-imaging ${D}${bindir}
		 cp -r ${WORKDIR}/data/ ${D}${bindir}
}
