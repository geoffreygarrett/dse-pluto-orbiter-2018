
# # DECORATORS ---------------------------------------------------------------------------------------------------------
# def decorators(handle_traffic=True, update_signature=False, retry_if_error=True):
#     class _Decorator(object):
#
#         def __init__(self, decorated):
#             self._decorated = decorated
#
#         def __get__(self, instance, owner):
#             self.cls = owner
#             self.obj = instance
#             return self.__call__
#
#         # Decorator 1
#         def _update_signature(self):
#
#         # Decorator 2
#         def _handle_traffic(self):
#
#         def __call__(self, *args, **kwargs):
#             if handle_traffic is True:
#                 self._handle_traffic()
#
#             if update_signature is True:
#                 self._update_signature()
#
#     return _Decorator