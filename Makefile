include $(TOPDIR)/rules.mk

PKG_NAME:=rtlsdr_wsprd
PKG_VERSION:=0.2
PKG_RELEASE:=2.00

PKG_SOURCE_PROTO:=git
PKG_SOURCE_URL:=https://github.com/Guenael/rtlsdr-wsprd.git
PKG_SOURCE_SUBDIR:=$(PKG_NAME)-$(PKG_VERSION)
PKG_SOURCE_VERSION:=master
PKG_SOURCE:=$(PKG_NAME)-$(PKG_VERSION)-$(PKG_SOURCE_VERSION).tar.gz

PKG_MAINTAINER:=Guenael
#PKG_BUILD_DIR:=$(BUILD_DIR)/$(PKG_NAME)/$(BUILD_VARIANT)/$(PKG_NAME)-$(PKG_VERSION)
PKG_BUILD_DIR:=$(BUILD_DIR)/$(PKG_SOURCE_SUBDIR)

include $(INCLUDE_DIR)/package.mk

define Package/rtlsdr_wsprd/Default
  TITLE:=This non-interactive application allows automatic reporting of WSPR spots on WSPRnet. \
It is written in C++ and uses RTL-SDR to interface with RTL2832-based \
hardware.
  URL:=https://github.com/Guenael/rtlsdr-wsprd.git
endef

define Build/Prepare
	mkdir -p $(PKG_BUILD_DIR)
	$(CP) /home/openwrt/openwrt_rtlsdr-wsprd/src/* $(PKG_BUILD_DIR)/
#	$(call Build/Prepare/Default)
#	export PKG_CONFIG_PATH=$(BUILD_DIR)/rtl-sdr-0.5.3-20150608-/src/ 
endef

#TARGET_LDFLAGS+= -L$(TOOLCHAIN_DIR)/usr/lib -L$(TOOLCHAIN_DIR)/lib -Wl,-rpath=$(TOOLCHAIN_DIR)/lib

define Package/rtlsdr_wsprd
  $(call Package/rtlsdr_wsprd/Default)
  SECTION:=utils
  CATEGORY:=Utilities
  DEPENDS:=+librtlsdr +libcurl +ntpdate +libpthread +fftw3
endef

define Package/rtlsdr_wsprd/description
  The idea is to allow the use of small computer like RasberryPi or Beaglebone boards, \
  with a simple deamon. This kind of very lightweight setup could run continuously \
  without maintenance and help to increase the WSPR network
endef

#define Package/csdr/conffiles
#  /etc/config/csdr
#endef

define Package/rtlsdr_wsprd/install
	$(INSTALL_DIR) $(1)/usr/bin
	$(INSTALL_BIN) $(PKG_BUILD_DIR)/$(PKG_NAME) $(1)/usr/bin 
#	$(INSTALL_DIR) $(1)/etc/multimon-ng
#	$(INSTALL_DATA) ./files/multimon-ng.template $(1)/etc/multimon-ng/config.template
endef

$(eval $(call BuildPackage,rtlsdr_wsprd))

