from django.conf.urls.defaults import *

# Uncomment the next two lines to enable the admin:
#from django.contrib import admin
#admin.autodiscover()

urlpatterns = patterns('',
    # Example:
    # (r'^LeafWebApp/', include('LeafWebApp.foo.urls')),
    (r'^', include('Leaf.urls')),
#    (r'^TEST/$', 'Leaf.TEST.index'),
#    (r'^TEST2/$', 'Leaf.TEST2.index'),

    # Uncomment the admin/doc line below and add 'django.contrib.admindocs' 
    # to INSTALLED_APPS to enable admin documentation:
    # (r'^admin/doc/', include('django.contrib.admindocs.urls')),

    # Uncomment the next line to enable the admin:
#    (r'^admin/', include(admin.site.urls)),
)
