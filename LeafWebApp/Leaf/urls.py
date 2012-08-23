from django.conf.urls.defaults import patterns, include, url

urlpatterns = patterns('',

    url(r'^$', 'Leaf.views.index'),
    url(r'^contact/$', 'Leaf.views.contact'),
    url(r'^culling/$', 'Leaf.views.culling'),
    url(r'^cullingchoice/$', 'Leaf.views.cullingchoice'),
    url(r'^user_culling/$', 'Leaf.views.user_culling'),
    url(r'^user_pdb_culling/$', 'Leaf.views.user_pdb_culling'),
    url(r'^whole_pdb_culling/$', 'Leaf.views.whole_pdb_culling'),
    url(r'^helppage/$', 'Leaf.views.help_page'),
    url(r'^culling/user_submit/$', 'Leaf.views.user_submit'),
    url(r'^culling/user_pdb_submit/$', 'Leaf.views.user_pdb_submit'),
#    url(r'^culling/whole_pdb_submit/$', 'Leaf.views.whole_pdb_submit'),
    url(r'^downloads/$', 'Leaf.views.downloads'),
    url(r'^downloads/(?P<fileName>[0-9a-zA-Z_./]+)/$', 'Leaf.views.download_gzipped'),
#    url(r'^requestsent/$', 'Leaf.views.sent'),
    url(r'^requestsent/(?P<result_id>\d+)/(?P<cullType>[a-z]+)/$', 'Leaf.views.sent'),
#    url(r'^results/[0-9]+/[0-9]+/[0-9]+/(?P<result_id>\d+)/(?P<cullType>[a-z]+)/$', 'Leaf.views.results'),
    url(r'^results/[0-9]+/[0-9]+/[0-9]+/(?P<result_id>\d+)/(?P<cullType>[a-z]+)/(?P<fileName>[a-zA-Z_]+)/$', 'Leaf.views.txtlist')
)
